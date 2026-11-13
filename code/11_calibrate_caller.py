#!/usr/bin/env python3
"""Calibrate the signal-to-peak caller against published MACS2 peak calls.

GSE162976 published, for the same GSMs, both a normalised signal track and an
ENCODE-pipeline MACS2 peak set. That pairing is used here as ground truth: we
sweep (quantile, gap, min_width), score each setting against the published peaks
at both base-pair and peak level, and record the setting that maximises agreement.
The chosen setting is then applied unchanged to the two cohorts that released
signal only, so their peak calls are not tuned on themselves.

Outputs
  results/calibration_grid.csv   every sample x parameter combination
  results/calibration_best.json  the selected parameters and their scores
"""
import importlib.util
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REF = os.path.join(DATA, "ref")
RES = os.path.join(ROOT, "results")
CACHE = os.path.join(DATA, "cache_bins")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sp = load_module("signal_peaks", os.path.join(HERE, "10_signal_peaks.py"))

QUANTILES = [0.970, 0.980, 0.985, 0.990, 0.995]
GAPS = [200, 400, 800]
MIN_WIDTHS = [500, 1000]


def mask_from_intervals(intervals, sizes, binsize=sp.BIN):
    """Boolean per-bin mask per chromosome."""
    out = {c: np.zeros(max(1, n // binsize), dtype=bool) for c, n in sizes.items()}
    for c, ivs in intervals.items():
        if c not in out:
            continue
        m = out[c]
        for s, e in ivs:
            a, b = s // binsize, min(len(m), (e + binsize - 1) // binsize)
            if b > a:
                m[a:b] = True
    return out


def jaccard(m1, m2):
    inter = union = 0
    for c in m1:
        if c not in m2:
            continue
        a, b = m1[c], m2[c]
        n = min(len(a), len(b))
        inter += int(np.count_nonzero(a[:n] & b[:n]))
        union += int(np.count_nonzero(a[:n] | b[:n]))
    return inter / union if union else 0.0


def peak_level(called, truth_iv, sizes):
    """precision = called regions hitting a published peak; recall = reverse."""
    tmask = mask_from_intervals(truth_iv, sizes)
    hit = 0
    for c, s, e, _, _ in called:
        m = tmask.get(c)
        if m is None:
            continue
        a, b = s // sp.BIN, min(len(m), (e + sp.BIN - 1) // sp.BIN)
        if b > a and m[a:b].any():
            hit += 1
    precision = hit / len(called) if called else 0.0

    cmask = mask_from_intervals(
        {c: [(s, e)] for c, s, e, _, _ in called} if False else _group(called), sizes
    )
    tot = cov = 0
    for c, ivs in truth_iv.items():
        m = cmask.get(c)
        for s, e in ivs:
            tot += 1
            if m is None:
                continue
            a, b = s // sp.BIN, min(len(m), (e + sp.BIN - 1) // sp.BIN)
            if b > a and m[a:b].any():
                cov += 1
    recall = cov / tot if tot else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return precision, recall, f1


def _group(called):
    d = {}
    for c, s, e, _, _ in called:
        d.setdefault(c, []).append((s, e))
    return d


def main():
    os.makedirs(RES, exist_ok=True)
    sizes = sp.read_sizes(os.path.join(REF, "hg38.main.chrom.sizes"))
    blacklist = os.path.join(REF, "hg38.blacklist.bed")

    mf = pd.read_csv(os.path.join(DATA, "manifest.csv"))
    ped = mf[(mf.cohort == "pHGG_pediatric") & (mf.mark == "H3K27ac")]
    sig = ped[(ped["kind"] == "signal") & ped["present"]].set_index("gsm")
    pk = ped[(ped["kind"] == "peaks") & ped["present"]].set_index("gsm")
    shared = sorted(set(sig.index) & set(pk.index))
    print(f"calibration samples with signal+peaks on disk: {len(shared)}")
    if not shared:
        print("no paired samples downloaded yet - rerun when the pedi bigwigs land")
        return 1

    rows = []
    for i, gsm in enumerate(shared, 1):
        sbw, pbed = sig.loc[gsm, "path"], pk.loc[gsm, "path"]
        label = sig.loc[gsm, "label"]
        truth = sp.load_bed_intervals(pbed)
        ntruth = sum(len(v) for v in truth.values())
        binned = sp.bin_genome(sbw, sizes, cache_dir=CACHE)
        pool = sp.nonzero_pool(binned)
        tmask = mask_from_intervals(truth, sizes)
        print(f"[{i}/{len(shared)}] {label}: {ntruth} published peaks", flush=True)
        for q, gap, mw in itertools.product(QUANTILES, GAPS, MIN_WIDTHS):
            called = sp.call_from_binned(binned, q, gap, mw, pool=pool)
            called = sp.drop_blacklist(called, blacklist)
            cmask = mask_from_intervals(_group(called), sizes)
            jac = jaccard(cmask, tmask)
            prec, rec, f1 = peak_level(called, truth, sizes)
            rows.append(
                {
                    "gsm": gsm, "label": label, "quantile": q, "gap": gap,
                    "min_width": mw, "n_called": len(called), "n_truth": ntruth,
                    "mb_called": sum(e - s for _, s, e, _, _ in called) / 1e6,
                    "jaccard": jac, "precision": prec, "recall": rec, "f1": f1,
                }
            )
        del binned

    grid = pd.DataFrame(rows)
    grid.to_csv(os.path.join(RES, "calibration_grid.csv"), index=False)

    agg = (
        grid.groupby(["quantile", "gap", "min_width"])
        .agg(
            jaccard_med=("jaccard", "median"),
            jaccard_iqr=("jaccard", lambda s: s.quantile(0.75) - s.quantile(0.25)),
            f1_med=("f1", "median"),
            precision_med=("precision", "median"),
            recall_med=("recall", "median"),
            n_called_med=("n_called", "median"),
            mb_med=("mb_called", "median"),
        )
        .reset_index()
        .sort_values("jaccard_med", ascending=False)
    )
    with pd.option_context("display.width", 200, "display.max_rows", 40):
        print("\n== top 12 parameter settings (median over samples) ==")
        print(agg.head(12).round(4).to_string(index=False))

    best = agg.iloc[0]
    sel = {
        "quantile": float(best["quantile"]),
        "gap": int(best["gap"]),
        "min_width": int(best["min_width"]),
        "binsize": sp.BIN,
        "jaccard_median": float(best["jaccard_med"]),
        "f1_median": float(best["f1_med"]),
        "precision_median": float(best["precision_med"]),
        "recall_median": float(best["recall_med"]),
        "n_called_median": float(best["n_called_med"]),
        "n_calibration_samples": len(shared),
        "calibration_cohort": "pHGG_pediatric (GSE162976)",
        "truth": "published ENCODE-pipeline MACS2 narrowPeak, same GSMs",
    }
    with open(os.path.join(RES, "calibration_best.json"), "w") as fh:
        json.dump(sel, fh, indent=2)
    print("\nselected:", json.dumps(sel, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
