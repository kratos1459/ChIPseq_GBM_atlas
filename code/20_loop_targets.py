#!/usr/bin/env python3
"""Enhancer-to-gene assignment: 3D contacts versus nearest-TSS.

Almost every published enhancer catalogue assigns each element to its nearest
gene. GSE262089 provides H3K27ac HiChIP loops (FitHiChIP, 10 kb anchors) for
four glioblastoma stem-cell lines and a neural stem-cell control, so the
nearest-gene assumption can be measured rather than assumed.

For each sample:
  * classify loop anchors as promoter (contains a TSS), enhancer (contains an
    H3K27ac peak but no TSS), or neither;
  * take every enhancer-promoter loop as a contact-based link;
  * for each linked enhancer, ask whether its nearest TSS is the gene it
    actually contacts, how far the true target is, and how many TSSs sit
    between the enhancer and its contacted gene.

Outputs (results/loop_targets/)
  links_<sample>.tsv          enhancer-promoter links with both assignments
  anchor_classes.csv          anchor composition per sample
  nearest_vs_contact.csv      per-sample agreement statistics
  hub_genes.csv               genes contacted by many distinct enhancers
"""
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
REF = os.path.join(DATA, "ref")
OUT = os.path.join(ROOT, "results", "loop_targets")

# GSE262089: loop file stem -> peak file stem (lab used different prefixes)
PAIRS = {
    "s320": "gsc320", "s728": "gsc728", "s810": "gsc810",
    "s1206": "gsc1206", "nsc": "nsc",
}


def read_bed3(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            if len(f) >= 3 and "_" not in f[0] and f[0] != "chrM":
                rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def read_bedpe(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0].startswith("#"):
                continue
            try:
                rows.append(
                    (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]),
                     float(f[6]) if len(f) > 6 and f[6] not in (".", "") else np.nan)
                )
            except ValueError:
                continue
    return rows


def read_tss(path):
    by_chrom = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or "_" in f[0]:
                continue
            by_chrom[f[0]].append((int(f[1]), f[3], f[6]))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def index_intervals(rows):
    d = defaultdict(list)
    for c, s, e in rows:
        d[c].append((s, e))
    for c in d:
        d[c].sort()
    return d


def overlaps(index, chrom, start, end):
    ivs = index.get(chrom)
    if not ivs:
        return False
    lo, hi = 0, len(ivs)
    while lo < hi:                      # first interval with end > start
        mid = (lo + hi) // 2
        if ivs[mid][1] <= start:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(ivs) and ivs[lo][0] < end


def tss_in(tss, chrom, start, end):
    arr = tss.get(chrom)
    if not arr:
        return []
    pos = [a[0] for a in arr]
    lo = np.searchsorted(pos, start, "left")
    hi = np.searchsorted(pos, end, "right")
    return arr[lo:hi]


def nearest_tss(tss, chrom, point):
    arr = tss.get(chrom)
    if not arr:
        return None, None, None
    pos = np.array([a[0] for a in arr])
    i = int(np.searchsorted(pos, point))
    best, bestd = None, None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(arr):
            d = abs(arr[j][0] - point)
            if bestd is None or d < bestd:
                best, bestd = arr[j], d
    return (best[1], bestd, best[0]) if best else (None, None, None)


def tss_between(tss, chrom, a, b):
    """How many distinct genes have a TSS strictly between two points."""
    lo, hi = min(a, b), max(a, b)
    return len({g for _, g, _ in tss_in(tss, chrom, lo, hi)})


def main():
    os.makedirs(OUT, exist_ok=True)
    tss = read_tss(os.path.join(REF, "gencode.v44.tss.hg38.bed"))
    mf = pd.read_csv(os.path.join(DATA, "manifest.csv"))
    g = mf[(mf.gse == "GSE262089") & mf.present]
    loops = {}
    peaks = {}
    for _, r in g.iterrows():
        f = r["file"]
        if f.endswith(".bedpe.gz"):
            stem = f.split("_")[1].replace(".bedpe.gz", "")
            loops[stem] = r["path"]
        elif f.endswith(".bed.gz") and "atac_" not in f:
            stem = f.split("_")[1].replace(".bed.gz", "")
            peaks[stem] = r["path"]

    print("loop files:", sorted(loops))
    print("peak files:", sorted(peaks))

    anchor_rows, stat_rows, hub_rows = [], [], []
    for lstem, lpath in sorted(loops.items()):
        pstem = PAIRS.get(lstem)
        ppath = peaks.get(pstem)
        if ppath is None:
            print(f"  ! no H3K27ac peaks for {lstem}; skipped")
            continue
        sample = lstem
        pk = index_intervals(read_bed3(ppath))
        lp = read_bedpe(lpath)
        npk = sum(len(v) for v in pk.values())
        print(f"\n== {sample}: {len(lp)} loops, {npk} H3K27ac peaks", flush=True)

        def classify(chrom, s, e):
            has_tss = bool(tss_in(tss, chrom, s, e))
            has_pk = overlaps(pk, chrom, s, e)
            if has_tss:
                return "promoter"
            return "enhancer" if has_pk else "other"

        cls_count = defaultdict(int)
        links = []
        for c1, s1, e1, c2, s2, e2, score in lp:
            k1, k2 = classify(c1, s1, e1), classify(c2, s2, e2)
            cls_count[(k1, k2) if k1 <= k2 else (k2, k1)] += 1
            pairs = []
            if k1 == "enhancer" and k2 == "promoter":
                pairs.append(((c1, s1, e1), (c2, s2, e2)))
            elif k2 == "enhancer" and k1 == "promoter":
                pairs.append(((c2, s2, e2), (c1, s1, e1)))
            for (ec, es, ee), (pc, ps, pe) in pairs:
                emid = (es + ee) // 2
                genes = sorted({gname for _, gname, _ in tss_in(tss, pc, ps, pe)})
                n_prom_tss = len(genes)
                nn, nd, npos = nearest_tss(tss, ec, emid)
                for gname in genes:
                    gpos = [p for p, gn, _ in tss_in(tss, pc, ps, pe) if gn == gname]
                    gp = gpos[0] if gpos else (ps + pe) // 2
                    links.append(
                        {
                            "sample": sample,
                            "enh_chrom": ec, "enh_start": es, "enh_end": ee,
                            "contact_gene": gname,
                            "contact_dist": abs(gp - emid),
                            "nearest_gene": nn,
                            "nearest_dist": nd,
                            "agree": int(nn == gname),
                            "genes_skipped": tss_between(tss, ec, emid, gp),
                            "n_prom_tss": n_prom_tss,
                            "loop_score": score,
                        }
                    )
        for (a, b), n in sorted(cls_count.items(), key=lambda x: -x[1]):
            anchor_rows.append({"sample": sample, "anchor_pair": f"{a}-{b}", "n_loops": n})

        if not links:
            print("  no enhancer-promoter loops found")
            continue
        df = pd.DataFrame(links)
        df.to_csv(os.path.join(OUT, f"links_{sample}.tsv"), sep="\t", index=False)

        enh = df.groupby(["enh_chrom", "enh_start", "enh_end"])
        per_enh_agree = enh["agree"].max()      # nearest gene is among contacts
        stat_rows.append(
            {
                "sample": sample,
                "n_loops": len(lp),
                "n_ep_links": len(df),
                "n_enhancers_linked": enh.ngroups,
                "n_genes_linked": df["contact_gene"].nunique(),
                "frac_nearest_is_a_contact": float(per_enh_agree.mean()),
                "frac_links_nearest_correct": float(df["agree"].mean()),
                "median_contact_dist_kb": float(df["contact_dist"].median() / 1e3),
                "median_nearest_dist_kb": float(df["nearest_dist"].median() / 1e3),
                "median_genes_skipped": float(df["genes_skipped"].median()),
                "mean_genes_skipped": float(df["genes_skipped"].mean()),
                "frac_links_skip_ge1_gene": float((df["genes_skipped"] >= 1).mean()),
                "enhancers_per_gene": float(len(df) / max(1, df["contact_gene"].nunique())),
            }
        )
        print(
            f"  E-P links {len(df)}, enhancers {enh.ngroups}, genes "
            f"{df['contact_gene'].nunique()}, nearest-gene correct "
            f"{df['agree'].mean():.1%}"
        )

        hub = (
            df.groupby("contact_gene")
            .agg(n_enhancers=("enh_start", "nunique"),
                 median_dist_kb=("contact_dist", lambda s: s.median() / 1e3))
            .sort_values("n_enhancers", ascending=False)
            .reset_index()
        )
        hub["sample"] = sample
        hub_rows.append(hub)

    if stat_rows:
        st = pd.DataFrame(stat_rows)
        st.to_csv(os.path.join(OUT, "nearest_vs_contact.csv"), index=False)
        pd.DataFrame(anchor_rows).to_csv(os.path.join(OUT, "anchor_classes.csv"), index=False)
        hubs = pd.concat(hub_rows, ignore_index=True)
        hubs.to_csv(os.path.join(OUT, "hub_genes.csv"), index=False)
        with pd.option_context("display.width", 220):
            print("\n== nearest-TSS versus 3D contact ==")
            print(st.round(3).to_string(index=False))
            print("\n== genes with most distinct contacting enhancers (GSC samples) ==")
            gsc = hubs[hubs["sample"] != "nsc"]
            top = (
                gsc.groupby("contact_gene")
                .agg(samples=("sample", "nunique"), enh=("n_enhancers", "sum"))
                .sort_values(["samples", "enh"], ascending=False)
                .head(25)
            )
            print(top.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
