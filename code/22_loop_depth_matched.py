#!/usr/bin/env python3
"""Depth-matched version of the GSC-versus-NSC contact comparison.

The five HiChIP samples differ several-fold in loop count (33k for NSC against
96k for GSC 320), so any gene that looks "contacted in every GSC line but not in
NSC" may simply reflect deeper loop calling. Here every sample is subsampled to
the same number of enhancer-promoter loops as the shallowest sample, repeatedly,
and the comparison is recomputed.

Outputs (results/loop_targets/)
  depth_matched_summary.csv     statistics at matched depth, mean over repeats
  depth_matched_genes.csv       genes still GSC-specific in most repeats
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "loop_targets")
N_REPEAT = 20
RNG = np.random.default_rng(20261113)


def loop_key(df):
    """One row per (enhancer anchor, promoter anchor) loop."""
    return list(zip(df["enh_chrom"], df["enh_start"], df["prom_start"]))


def main():
    files = sorted(glob.glob(os.path.join(OUT, "links_*.tsv")))
    if not files:
        print("no links_*.tsv - run 20_loop_targets.py first")
        return 1
    frames = {}
    for f in files:
        s = os.path.basename(f)[len("links_"):-len(".tsv")]
        d = pd.read_csv(f, sep="\t")
        if "prom_start" not in d.columns:
            print("links files predate the prom_start column - re-run "
                  "20_loop_targets.py so loops can be identified")
            return 1
        d["_loop"] = loop_key(d)
        frames[s] = d
    gsc = sorted(s for s in frames if s != "nsc")
    if "nsc" not in frames or not gsc:
        print("need both GSC and NSC samples")
        return 1

    loops = {s: pd.unique(pd.Series(d["_loop"])) for s, d in frames.items()}
    counts = {s: len(v) for s, v in loops.items()}
    print("distinct enhancer-promoter loops:", counts)
    target = min(counts.values())
    print(f"matching every sample to {target} loops, {N_REPEAT} repeats")

    rows, gene_hits = [], {}
    for rep in range(N_REPEAT):
        kept = {}
        for s, d in frames.items():
            sel = set(map(tuple, np.asarray(
                loops[s], dtype=object)[RNG.choice(counts[s], target, replace=False)]
            )) if counts[s] > target else set(map(tuple, loops[s]))
            kept[s] = d[[tuple(k) in sel for k in d["_loop"]]]
        gsets = {s: set(k["contact_gene"]) for s, k in kept.items()}
        spec = set.intersection(*[gsets[s] for s in gsc]) - gsets["nsc"]
        nsc_spec = gsets["nsc"] - set.union(*[gsets[s] for s in gsc])
        for g in spec:
            gene_hits[g] = gene_hits.get(g, 0) + 1
        rows.append({
            "repeat": rep,
            "genes_all_gsc_not_nsc": len(spec),
            "genes_nsc_not_any_gsc": len(nsc_spec),
            "genes_contacted_nsc": len(gsets["nsc"]),
            "genes_contacted_gsc_mean": float(np.mean([len(gsets[s]) for s in gsc])),
            "nearest_correct_gsc_mean": float(np.mean(
                [kept[s]["agree"].mean() for s in gsc])),
            "nearest_correct_nsc": float(kept["nsc"]["agree"].mean()),
        })
        if rep == 0:
            print(f"  repeat 0: {len(spec)} GSC-specific, {len(nsc_spec)} NSC-specific")

    summ = pd.DataFrame(rows)
    agg = summ.drop(columns=["repeat"]).agg(["mean", "std"]).T
    agg.columns = ["mean", "sd"]
    agg.to_csv(os.path.join(OUT, "depth_matched_summary.csv"))
    print("\n== at matched loop depth (mean +/- sd over repeats) ==")
    print(agg.round(3).to_string())

    gh = (
        pd.DataFrame({"gene": list(gene_hits), "n_repeats_specific": list(gene_hits.values())})
        .sort_values("n_repeats_specific", ascending=False)
    )
    gh["frac_repeats"] = gh["n_repeats_specific"] / N_REPEAT
    gh.to_csv(os.path.join(OUT, "depth_matched_genes.csv"), index=False)
    robust = gh[gh["frac_repeats"] >= 0.8]
    named = [g for g in robust["gene"] if not str(g).startswith("ENSG")]
    print(f"\ngenes GSC-specific in >=80% of depth-matched repeats: {len(robust)} "
          f"({len(named)} with symbols)")
    print("examples:", sorted(named)[:40])

    # how much of the unmatched claim survives matching
    before = os.path.join(OUT, "gsc_specific_targets.csv")
    if os.path.exists(before):
        b = pd.read_csv(before)
        n_before = int(((b["n_gsc_samples"] == len(gsc)) & (b["in_nsc"] == 0)).sum())
        print(f"\nunmatched analysis called {n_before} GSC-specific genes; "
              f"{len(robust)} survive depth matching "
              f"({len(robust)/n_before:.1%})" if n_before else "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
