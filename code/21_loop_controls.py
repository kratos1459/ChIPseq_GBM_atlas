#!/usr/bin/env python3
"""Controls and reproducibility for the contact-versus-nearest-gene result.

A low nearest-gene hit rate is only interesting if it is not an artefact of
(i) promoter anchors that contain several TSSs, or (ii) the trivial fact that a
gene 100 kb away is rarely the closest one. Both are measured here, together
with how reproducible contact-based targets are across independent GSC lines.

  1. clean subset      links whose promoter anchor contains exactly one TSS
  2. distance-matched  for every link, a decoy gene at the same distance from the
     null             same enhancer; how often is the decoy the nearest gene?
  3. reproducibility   Jaccard of contacted-gene sets for enhancers shared
                       between GSC lines, against a within-cohort noise ceiling
  4. GSC vs NSC        genes contacted in every GSC line but not in NSC

Outputs (results/loop_targets/)
  controls_summary.csv, reproducibility_pairs.csv, gsc_specific_targets.csv
"""
import glob
import importlib.util
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REF = os.path.join(DATA, "ref")
OUT = os.path.join(ROOT, "results", "loop_targets")
RNG = np.random.default_rng(20261113)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lt = load_module("loop_targets", os.path.join(HERE, "20_loop_targets.py"))


def decoy_agreement(df, tss, n_max=40000):
    """For each link, pick a gene at a matched distance from the same enhancer
    (excluding the true contact) and ask how often it is the nearest gene."""
    sub = df if len(df) <= n_max else df.sample(n_max, random_state=7)
    hits = tot = 0
    pos_cache = {}
    for r in sub.itertuples(index=False):
        c = r.enh_chrom
        if c not in pos_cache:
            arr = tss.get(c, [])
            pos_cache[c] = (
                np.array([a[0] for a in arr]),
                np.array([a[1] for a in arr], dtype=object),
            )
        pos, names = pos_cache[c]
        if pos.size < 3:
            continue
        emid = (r.enh_start + r.enh_end) // 2
        d = np.abs(pos - emid)
        order = np.argsort(np.abs(d - r.contact_dist))
        for j in order[:4]:
            if names[j] != r.contact_gene:
                hits += int(names[j] == r.nearest_gene)
                tot += 1
                break
    return hits / tot if tot else np.nan, tot


def main():
    files = sorted(glob.glob(os.path.join(OUT, "links_*.tsv")))
    if not files:
        print("no links_*.tsv - run 20_loop_targets.py first")
        return 1
    tss = lt.read_tss(os.path.join(REF, "gencode.v44.tss.hg38.bed"))

    frames = {}
    for f in files:
        s = os.path.basename(f)[len("links_"):-len(".tsv")]
        frames[s] = pd.read_csv(f, sep="\t")
    print("samples:", sorted(frames))

    rows = []
    for s, df in frames.items():
        clean = df[df["n_prom_tss"] == 1]
        null_rate, null_n = decoy_agreement(df, tss)
        enh = df.groupby(["enh_chrom", "enh_start", "enh_end"])
        rows.append(
            {
                "sample": s,
                "n_links": len(df),
                "nearest_correct_all": float(df["agree"].mean()),
                "n_links_clean": len(clean),
                "nearest_correct_clean": float(clean["agree"].mean()) if len(clean) else np.nan,
                "nearest_correct_decoy_null": null_rate,
                "decoy_n": null_n,
                "enh_with_nearest_among_contacts": float(enh["agree"].max().mean()),
                "median_contact_kb": float(df["contact_dist"].median() / 1e3),
                "median_nearest_kb": float(df["nearest_dist"].median() / 1e3),
                "targets_per_enhancer": float(enh["contact_gene"].nunique().mean()),
            }
        )
    summ = pd.DataFrame(rows).sort_values("sample")
    summ.to_csv(os.path.join(OUT, "controls_summary.csv"), index=False)
    with pd.option_context("display.width", 230):
        print("\n== controls ==")
        print(summ.round(4).to_string(index=False))

    # ---- reproducibility of contact targets between GSC lines ---------------
    def target_sets(df):
        d = {}
        for r in df.itertuples(index=False):
            d.setdefault((r.enh_chrom, r.enh_start), set()).add(r.contact_gene)
        return d

    sets = {s: target_sets(df) for s, df in frames.items()}
    gsc = [s for s in sets if s != "nsc"]
    pair_rows = []
    for a, b in itertools.combinations(sorted(gsc), 2):
        shared = set(sets[a]) & set(sets[b])
        if not shared:
            continue
        jac, same_any = [], []
        for k in shared:
            A, B = sets[a][k], sets[b][k]
            jac.append(len(A & B) / len(A | B))
            same_any.append(int(bool(A & B)))
        pair_rows.append(
            {
                "sample_a": a, "sample_b": b,
                "shared_enhancers": len(shared),
                "enh_a": len(sets[a]), "enh_b": len(sets[b]),
                "enhancer_overlap_frac": len(shared) / min(len(sets[a]), len(sets[b])),
                "mean_target_jaccard": float(np.mean(jac)),
                "median_target_jaccard": float(np.median(jac)),
                "frac_any_shared_target": float(np.mean(same_any)),
            }
        )
    # NSC comparison for contrast
    for a in sorted(gsc):
        shared = set(sets[a]) & set(sets["nsc"])
        if not shared:
            continue
        jac = [
            len(sets[a][k] & sets["nsc"][k]) / len(sets[a][k] | sets["nsc"][k])
            for k in shared
        ]
        pair_rows.append(
            {
                "sample_a": a, "sample_b": "nsc",
                "shared_enhancers": len(shared),
                "enh_a": len(sets[a]), "enh_b": len(sets["nsc"]),
                "enhancer_overlap_frac": len(shared) / min(len(sets[a]), len(sets["nsc"])),
                "mean_target_jaccard": float(np.mean(jac)),
                "median_target_jaccard": float(np.median(jac)),
                "frac_any_shared_target": float(np.mean([j > 0 for j in jac])),
            }
        )
    rep = pd.DataFrame(pair_rows)
    rep.to_csv(os.path.join(OUT, "reproducibility_pairs.csv"), index=False)
    with pd.option_context("display.width", 230):
        print("\n== target-set reproducibility between samples ==")
        print(rep.round(4).to_string(index=False))

    # ---- GSC-specific target genes -----------------------------------------
    gene_samples = {}
    for s, df in frames.items():
        for g in df["contact_gene"].unique():
            gene_samples.setdefault(g, set()).add(s)
    recs = []
    for g, ss in gene_samples.items():
        ngsc = len(ss & set(gsc))
        recs.append({"gene": g, "n_gsc_samples": ngsc, "in_nsc": int("nsc" in ss)})
    gs = pd.DataFrame(recs)
    spec = gs[(gs["n_gsc_samples"] == len(gsc)) & (gs["in_nsc"] == 0)]
    gs.to_csv(os.path.join(OUT, "gsc_specific_targets.csv"), index=False)
    print(
        f"\ngenes contacted in all {len(gsc)} GSC lines but not NSC: {len(spec)}"
        f"  (of {len(gs)} genes contacted anywhere)"
    )
    named = [x for x in spec["gene"] if not x.startswith("ENSG")]
    print("examples:", sorted(named)[:40])
    return 0


if __name__ == "__main__":
    sys.exit(main())
