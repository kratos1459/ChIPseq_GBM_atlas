#!/usr/bin/env python3
"""How much of the glioblastoma super-enhancer landscape reproduces between
independent cohorts, and is cross-cohort agreement near the ceiling set by
sampling noise?

The genome is divided into fixed 5 kb bins (a union-merge of stitched enhancer
domains across 160+ samples collapses into megabase blocks, so bins are used
instead). Every sample contributes two binary tracks: bins covered by a
super-enhancer, and bins covered by any stitched enhancer. Cohort recurrence is
then the fraction of samples in that cohort carrying a super-enhancer in a bin.

Agreement between two cohorts is only interpretable against the agreement two
halves of the *same* cohort reach, so a split-half noise ceiling is computed for
each cohort and reported alongside every cross-cohort statistic.

Outputs (results/cross_cohort/)
  recurrence.csv.gz        per-bin recurrence in every cohort group
  agreement.csv            cross-cohort agreement with split-half ceilings
  core_domains.bed         consensus GBM super-enhancer domains
  core_domain_genes.csv    nearest-TSS and HiChIP-contact targets per domain
  summary.json
"""
import glob
import gzip
import importlib.util
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats as sstats
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REF = os.path.join(DATA, "ref")
RES = os.path.join(ROOT, "results")
SEDIR = os.path.join(RES, "se")
OUT = os.path.join(RES, "cross_cohort")
BINSIZE = 5000
RNG = np.random.default_rng(20261113)
N_SPLITS = 25
NSC_MAX = 0.25       # at most 2 of 12 NSC samples (v2; was 0.10)
OTHER_MAX = 0.25     # 'absent' from the other adult cohort (v2; was 0.10)
TOP_N = 500          # fixed-cutoff sensitivity analysis: top domains per sample


def lm(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lt = lm("loop_targets", os.path.join(HERE, "20_loop_targets.py"))


class BinIndex:
    """Flat index over fixed-size bins of the main chromosomes."""

    def __init__(self, sizes_path, binsize=BINSIZE):
        sizes = {}
        with open(sizes_path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    sizes[f[0]] = int(f[1])
        self.binsize = binsize
        self.chroms = sorted(sizes, key=lambda c: (len(c), c))
        self.nbins, self.offset, self.counts = 0, {}, {}
        for c in self.chroms:
            n = sizes[c] // binsize + 1
            self.offset[c] = self.nbins
            self.counts[c] = n
            self.nbins += n

    def mark(self, mask, chrom, start, end):
        off = self.offset.get(chrom)
        if off is None:
            return
        a = off + start // self.binsize
        b = off + min(self.counts[chrom] - 1, (end - 1) // self.binsize) + 1
        if b > a:
            mask[a:b] = True

    def bin_to_locus(self, i):
        for c in self.chroms:
            o, n = self.offset[c], self.counts[c]
            if o <= i < o + n:
                s = (i - o) * self.binsize
                return c, s, s + self.binsize
        raise IndexError(i)

    def loci_frame(self, idx):
        chrom = np.empty(len(idx), dtype=object)
        start = np.zeros(len(idx), dtype=np.int64)
        bounds = [(c, self.offset[c], self.offset[c] + self.counts[c])
                  for c in self.chroms]
        for c, lo, hi in bounds:
            m = (idx >= lo) & (idx < hi)
            if m.any():
                chrom[m] = c
                start[m] = (idx[m] - lo) * self.binsize
        return pd.DataFrame({"chrom": chrom, "start": start,
                             "end": start + self.binsize})


def read_sample_bed(path):
    """Return (super-enhancers, all stitched domains, domains ranked by signal).

    The per-sample super-enhancer count depends on the shape of each sample's
    ranked-signal curve, so it varies several-fold between samples. The third
    return value supports a fixed top-N definition as a sensitivity analysis:
    if the cross-cohort conclusions hold under both an adaptive and a fixed
    cutoff, they are not an artefact of the cutoff rule.
    """
    se, enh, scored = [], [], []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 7 or "_" in f[0]:
                continue
            try:
                c, s, e, flag = f[0], int(f[1]), int(f[2]), int(f[6])
            except ValueError:
                continue
            sig = 0.0
            if "|" in f[3]:
                try:
                    sig = float(f[3].split("|", 1)[1])
                except ValueError:
                    sig = 0.0
            enh.append((c, s, e))
            scored.append((c, s, e, sig))
            if flag:
                se.append((c, s, e))
    scored.sort(key=lambda r: -r[3])
    return se, enh, scored


def group_of(cohort, tissue):
    if tissue == "normal/NSC":
        return "NSC"
    if cohort == "GSC_lines_adult":
        return "GSC_line"
    if cohort == "patient_tumour_adult":
        return "patient_tumour"
    if cohort == "pHGG_pediatric":
        return "pediatric"
    return "other"


def split_half_ceiling(mat, n_splits=N_SPLITS, thresh=0.5):
    """Agreement two random halves of one cohort reach on each other."""
    n = mat.shape[0]
    if n < 6:
        return {}
    out = {"spearman": [], "auc": [], "ap": [], "jaccard": []}
    for _ in range(n_splits):
        perm = RNG.permutation(n)
        a, b = perm[: n // 2], perm[n // 2 : 2 * (n // 2)]
        ra, rb = mat[a].mean(0), mat[b].mean(0)
        keep = (ra > 0) | (rb > 0)
        if keep.sum() < 100:
            continue
        out["spearman"].append(sstats.spearmanr(ra[keep], rb[keep]).statistic)
        yb = (rb[keep] >= thresh).astype(int)
        if 0 < yb.sum() < yb.size:
            out["auc"].append(roc_auc_score(yb, ra[keep]))
            out["ap"].append(average_precision_score(yb, ra[keep]))
        ya = (ra >= thresh)
        yb2 = (rb >= thresh)
        u = (ya | yb2).sum()
        out["jaccard"].append(((ya & yb2).sum() / u) if u else np.nan)
    return {k: (float(np.nanmean(v)), float(np.nanstd(v))) for k, v in out.items() if v}


def cross(ra, rb, thresh=0.5):
    keep = (ra > 0) | (rb > 0)
    res = {"n_bins_either": int(keep.sum())}
    if keep.sum() < 100:
        return res
    res["spearman"] = float(sstats.spearmanr(ra[keep], rb[keep]).statistic)
    yb = (rb[keep] >= thresh).astype(int)
    res["auc_predicting_b"] = (
        float(roc_auc_score(yb, ra[keep])) if 0 < yb.sum() < yb.size else np.nan
    )
    res["ap_predicting_b"] = (
        float(average_precision_score(yb, ra[keep])) if 0 < yb.sum() < yb.size else np.nan
    )
    res["ap_baseline"] = float(yb.mean())
    ya_t, yb_t = ra >= thresh, rb >= thresh
    u = (ya_t | yb_t).sum()
    res["jaccard_at_50pct"] = float((ya_t & yb_t).sum() / u) if u else np.nan
    res["n_a_recurrent"] = int(ya_t.sum())
    res["n_b_recurrent"] = int(yb_t.sum())
    res["frac_a_recurrent_also_b"] = (
        float((ya_t & yb_t).sum() / ya_t.sum()) if ya_t.sum() else np.nan
    )
    return res


def runs_to_bed(idx, bi):
    """Consecutive bin indices -> merged intervals."""
    if idx.size == 0:
        return []
    out = []
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    for s, e in zip(starts, ends):
        c1, s1, _ = bi.bin_to_locus(int(s))
        c2, _, e2 = bi.bin_to_locus(int(e))
        if c1 == c2:
            out.append((c1, s1, e2))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    bi = BinIndex(os.path.join(REF, "hg38.main.chrom.sizes"))
    print(f"{bi.nbins} bins of {BINSIZE} bp")

    mf = pd.read_csv(os.path.join(DATA, "manifest.csv"))
    meta = {}
    for _, r in mf[(mf["mark"] == "H3K27ac") & (mf["kind"] == "signal")].iterrows():
        lab = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(r["label"]))
        meta[(r["cohort"], lab)] = r

    files = sorted(glob.glob(os.path.join(SEDIR, "*", "*.hg38.bed")))
    print(f"{len(files)} per-sample region files")
    if len(files) < 10:
        print("too few samples processed yet - run 13_call_all.py first")
        return 1

    se_by_group, enh_by_group, labels = defaultdict(list), defaultdict(list), []
    top_by_group = defaultdict(list)
    for p in files:
        cohort = os.path.basename(os.path.dirname(p))
        lab = os.path.basename(p)[: -len(".hg38.bed")]
        r = meta.get((cohort, lab))
        tissue = r["tissue_class"] if r is not None else ""
        grp = group_of(cohort, tissue if isinstance(tissue, str) else "")
        se, enh, scored = read_sample_bed(p)
        m_se = np.zeros(bi.nbins, dtype=bool)
        m_en = np.zeros(bi.nbins, dtype=bool)
        m_top = np.zeros(bi.nbins, dtype=bool)
        for c, s, e in se:
            bi.mark(m_se, c, s, e)
        for c, s, e in enh:
            bi.mark(m_en, c, s, e)
        for c, s, e, _ in scored[:TOP_N]:
            bi.mark(m_top, c, s, e)
        se_by_group[grp].append(m_se)
        enh_by_group[grp].append(m_en)
        top_by_group[grp].append(m_top)
        labels.append({"cohort": cohort, "label": lab, "group": grp,
                       "se_bins": int(m_se.sum()), "enh_bins": int(m_en.sum()),
                       "n_se": len(se), "n_stitched": len(enh),
                       "top_bins": int(m_top.sum())})

    lab_df = pd.DataFrame(labels)
    lab_df.to_csv(os.path.join(OUT, "samples_used.csv"), index=False)
    print("\n== samples per group ==")
    print(lab_df.groupby("group").agg(n=("label", "size"),
                                      se_bins_med=("se_bins", "median")).to_string())

    mats = {g: np.vstack(v) for g, v in se_by_group.items() if len(v) >= 3}
    rec = {g: m.mean(0) for g, m in mats.items()}
    enh_mats = {g: np.vstack(v) for g, v in enh_by_group.items() if len(v) >= 3}
    rec_enh = {g: m.mean(0) for g, m in enh_mats.items()}

    rec_top_early = {g: np.vstack(v).mean(0) for g, v in top_by_group.items() if len(v) >= 3}
    any_se = np.zeros(bi.nbins, dtype=bool)
    for r in rec.values():
        any_se |= r > 0
    for r in rec_top_early.values():
        any_se |= r > 0
    idx = np.flatnonzero(any_se)
    loci = bi.loci_frame(idx)
    for g, r in rec.items():
        loci[f"se_rec_{g}"] = r[idx]
    for g, r in rec_enh.items():
        loci[f"enh_rec_{g}"] = r[idx]
    for g, r in rec_top_early.items():
        loci[f"top_rec_{g}"] = r[idx]
    loci.to_csv(os.path.join(OUT, "recurrence.csv.gz"), index=False,
                compression="gzip")
    print(f"\n{len(idx)} bins carry a super-enhancer in at least one sample "
          f"({len(idx)*BINSIZE/1e6:.0f} Mb)")

    # ---- ceilings and cross-cohort agreement -------------------------------
    rows = []
    ceil = {g: split_half_ceiling(m) for g, m in mats.items()}
    for g, c in ceil.items():
        if c:
            rows.append({
                "comparison": f"{g} split-half (ceiling)", "n_a": mats[g].shape[0],
                "n_b": mats[g].shape[0],
                "spearman": c.get("spearman", (np.nan,))[0],
                "spearman_sd": c.get("spearman", (np.nan, np.nan))[1],
                "auc_predicting_b": c.get("auc", (np.nan,))[0],
                "ap_predicting_b": c.get("ap", (np.nan,))[0],
                "auc_sd": c.get("auc", (np.nan, np.nan))[1],
                "jaccard_at_50pct": c.get("jaccard", (np.nan,))[0],
            })
    pairs = [("GSC_line", "patient_tumour"), ("patient_tumour", "GSC_line"),
             ("GSC_line", "pediatric"), ("patient_tumour", "pediatric"),
             ("GSC_line", "NSC"), ("patient_tumour", "NSC")]
    for a, b in pairs:
        if a in rec and b in rec:
            d = cross(rec[a], rec[b])
            d.update({"comparison": f"{a} -> {b}", "n_a": mats[a].shape[0],
                      "n_b": mats[b].shape[0]})
            rows.append(d)
    ag = pd.DataFrame(rows)
    ag["definition"] = "inflection (ROSE)"
    ag.to_csv(os.path.join(OUT, "agreement.csv"), index=False)
    with pd.option_context("display.width", 240):
        print("\n== agreement, with split-half ceilings ==")
        print(ag.round(4).to_string(index=False))

    # ---- sensitivity: fixed top-N cutoff instead of the inflection point ----
    top_mats = {g: np.vstack(v) for g, v in top_by_group.items() if len(v) >= 3}
    rec_top = {g: m.mean(0) for g, m in top_mats.items()}
    trows = []
    for g, m in top_mats.items():
        c = split_half_ceiling(m)
        if c:
            trows.append({
                "comparison": f"{g} split-half (ceiling)", "n_a": m.shape[0],
                "n_b": m.shape[0],
                "spearman": c.get("spearman", (np.nan,))[0],
                "auc_predicting_b": c.get("auc", (np.nan,))[0],
                "ap_predicting_b": c.get("ap", (np.nan,))[0],
                "jaccard_at_50pct": c.get("jaccard", (np.nan,))[0],
            })
    for a, b in pairs:
        if a in rec_top and b in rec_top:
            d = cross(rec_top[a], rec_top[b])
            d.update({"comparison": f"{a} -> {b}", "n_a": top_mats[a].shape[0],
                      "n_b": top_mats[b].shape[0]})
            trows.append(d)
    agt = pd.DataFrame(trows)
    agt["definition"] = f"fixed top-{TOP_N}"
    agt.to_csv(os.path.join(OUT, "agreement_topN.csv"), index=False)
    with pd.option_context("display.width", 240):
        print(f"\n== sensitivity: same statistics with a fixed top-{TOP_N} cutoff ==")
        print(agt.round(4).to_string(index=False))
    pd.concat([ag, agt], ignore_index=True).to_csv(
        os.path.join(OUT, "agreement_both_definitions.csv"), index=False)

    # ---- consensus core domains -------------------------------------------
    summary = {"binsize": BINSIZE, "n_bins_any_se": int(len(idx)),
               "definition": f"fixed top-{TOP_N} domains per sample",
               "nsc_max": NSC_MAX, "other_max": OTHER_MAX}
    R = rec_top
    if "GSC_line" in R and "patient_tumour" in R:
        nsc = R.get("NSC", np.zeros(bi.nbins))
        core = (R["GSC_line"] >= 0.5) & (R["patient_tumour"] >= 0.5)
        core_spec = core & (nsc < NSC_MAX)
        beds = runs_to_bed(np.flatnonzero(core_spec), bi)
        with open(os.path.join(OUT, "core_domains.bed"), "w") as fh:
            for i, (c, s, e) in enumerate(beds, 1):
                fh.write(f"{c}\t{s}\t{e}\tCORE_{i}\t0\t.\n")

        # contrast sets for the motif stage: cohort-specific and sporadic
        n_gsc, n_pat = mats["GSC_line"].shape[0], mats["patient_tumour"].shape[0]
        sets = {
            "gsc_specific": (R["GSC_line"] >= 0.5) & (R["patient_tumour"] < OTHER_MAX),
            "patient_specific": (R["patient_tumour"] >= 0.5) & (R["GSC_line"] < OTHER_MAX),
            "sporadic": (
                (R["GSC_line"] * n_gsc + R["patient_tumour"] * n_pat <= 1.001)
                & ((R["GSC_line"] > 0) | (R["patient_tumour"] > 0))
            ),
        }
        for name, mask in sets.items():
            b = runs_to_bed(np.flatnonzero(mask), bi)
            with open(os.path.join(OUT, f"{name}_domains.bed"), "w") as fh:
                for i, (c, s, e) in enumerate(b, 1):
                    fh.write(f"{c}\t{s}\t{e}\t{name.upper()}_{i}\t0\t.\n")
            summary[f"n_{name}_domains"] = len(b)
            summary[f"{name}_mb"] = round(int(mask.sum()) * BINSIZE / 1e6, 2)
            print(f"  {name}: {len(b)} domains, {mask.sum()*BINSIZE/1e6:.1f} Mb")
        summary.update({
            "n_bins_core_both_cohorts": int(core.sum()),
            "n_bins_core_tumour_specific": int(core_spec.sum()),
            "n_core_domains": len(beds),
            "core_mb": round(core_spec.sum() * BINSIZE / 1e6, 2),
        })
        print(f"\ncore domains (SE in >=50% of both adult cohorts, <25% of NSC): "
              f"{len(beds)} domains, {core_spec.sum()*BINSIZE/1e6:.1f} Mb")
        if "pediatric" in R:
            ped = R["pediatric"][core_spec]
            summary["core_median_recurrence_in_pediatric"] = float(np.median(ped))
            summary["core_frac_recurrent_in_pediatric"] = float((ped >= 0.5).mean())
            print(f"  of those, {(ped >= 0.5).mean():.1%} are also recurrent in "
                  f"pediatric pHGG (median recurrence {np.median(ped):.2f})")

        # ---- genes at core domains ----------------------------------------
        tss = lt.read_tss(os.path.join(REF, "gencode.v44.tss.hg38.bed"))
        links = {}
        for f in glob.glob(os.path.join(RES, "loop_targets", "links_*.tsv")):
            s = os.path.basename(f)[len("links_"):-len(".tsv")]
            if s == "nsc":
                continue
            d = pd.read_csv(f, sep="\t")
            for r in d.itertuples(index=False):
                links.setdefault((r.enh_chrom, r.enh_start // 10000), set()).add(
                    r.contact_gene
                )
        grows = []
        for i, (c, s, e) in enumerate(beds, 1):
            mid = (s + e) // 2
            nn, nd, _ = lt.nearest_tss(tss, c, mid)
            contacts = set()
            for b in range(s // 10000, e // 10000 + 1):
                contacts |= links.get((c, b), set())
            grows.append({
                "domain": f"CORE_{i}", "chrom": c, "start": s, "end": e,
                "width_kb": (e - s) / 1e3, "nearest_gene": nn,
                "nearest_dist_kb": (nd / 1e3) if nd is not None else np.nan,
                "n_contact_genes": len(contacts),
                "contact_genes": ";".join(sorted(x for x in contacts
                                                 if not x.startswith("ENSG"))[:12]),
                "nearest_is_contact": int(nn in contacts) if contacts else -1,
            })
        gd = pd.DataFrame(grows)
        gd.to_csv(os.path.join(OUT, "core_domain_genes.csv"), index=False)
        withc = gd[gd["nearest_is_contact"] >= 0]
        if len(withc):
            summary["core_nearest_is_contact_frac"] = float(
                withc["nearest_is_contact"].mean()
            )
            print(f"  core domains with loop data: {len(withc)}; nearest gene is a "
                  f"contact for {withc['nearest_is_contact'].mean():.1%}")
        print("\n== widest core domains ==")
        print(gd.nlargest(20, "width_kb")[
            ["domain", "chrom", "start", "width_kb", "nearest_gene", "contact_genes"]
        ].to_string(index=False))

    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\nsummary:", json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
