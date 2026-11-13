#!/usr/bin/env python3
"""Call enhancers and super-enhancers for every H3K27ac sample, one method for
all four cohorts.

Peak calling uses the parameters chosen by 11_calibrate_caller.py. Calling and
stitching happen in each cohort's native genome build; the resulting enhancer
coordinates are then lifted to hg38 so cohorts can be compared. Published peak
calls, where they exist, are left out of this step on purpose - every cohort is
processed identically so that cross-cohort differences cannot come from the
peak caller.

Resumable: a sample whose outputs already exist is skipped.

Outputs
  data/peaks_called/<cohort>/<label>.bed         enriched regions, native build
  results/se/<cohort>/<label>.tsv                stitched enhancers + SE flag
  results/se/<cohort>/<label>.hg38.bed           same regions lifted to hg38
  results/se_summary.csv                         one row per sample
"""
import argparse
import concurrent.futures as cf
import importlib.util
import json
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REF = os.path.join(DATA, "ref")
RES = os.path.join(ROOT, "results")
CACHE = os.path.join(DATA, "cache_bins")
PEAKS = os.path.join(DATA, "peaks_called")
SEDIR = os.path.join(RES, "se")
DEFAULTS = {"quantile": 0.99, "gap": 400, "min_width": 500}


def lm(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sp = lm("signal_peaks", os.path.join(HERE, "10_signal_peaks.py"))
se = lm("stitch_se", os.path.join(HERE, "12_stitch_superenhancers.py"))


def params():
    p = os.path.join(RES, "calibration_best.json")
    if os.path.exists(p):
        with open(p) as fh:
            d = json.load(fh)
        return {k: d[k] for k in ("quantile", "gap", "min_width")}, d.get(
            "jaccard_median"
        )
    print("! no calibration_best.json - falling back to defaults", file=sys.stderr)
    return dict(DEFAULTS), None


def safe(label):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label))


def signal_path(row):
    """Converted bigWig for the hg19 bedGraph cohort, raw file otherwise."""
    if row["gse"] == "GSE119776":
        base = os.path.basename(row["path"]).replace(".bedGraph.gz", ".bw")
        return os.path.join(DATA, "bigwig_hg19", base)
    return row["path"]


def liftover(in_bed, out_bed, chain):
    unmapped = out_bed + ".unmapped"
    p = subprocess.run(
        ["liftOver", in_bed, chain, out_bed, unmapped],
        capture_output=True, text=True,
    )
    return p.returncode == 0 and os.path.exists(out_bed)


def one_sample(row, par):
    cohort, label, genome = row["cohort"], safe(row["label"]), row["genome"]
    bw = signal_path(row)
    if not (bw and os.path.exists(bw)):
        return {"label": row["label"], "cohort": cohort, "status": "signal-missing"}

    peak_bed = os.path.join(PEAKS, cohort, f"{label}.bed")
    se_tsv = os.path.join(SEDIR, cohort, f"{label}.tsv")
    os.makedirs(os.path.dirname(peak_bed), exist_ok=True)
    os.makedirs(os.path.dirname(se_tsv), exist_ok=True)

    sizes = sp.read_sizes(os.path.join(REF, f"{genome}.main.chrom.sizes"))
    blacklist = os.path.join(REF, f"{genome}.blacklist.bed")
    promoters = os.path.join(REF, f"promoters_2kb.{genome}.bed")

    binned = None
    try:
        if not (os.path.exists(peak_bed) and os.path.exists(se_tsv)):
            binned = sp.bin_genome(bw, sizes, cache_dir=None)
        if not os.path.exists(peak_bed):
            regs = sp.call_from_binned(
                binned, par["quantile"], par["gap"], par["min_width"]
            )
            regs = sp.drop_blacklist(regs, blacklist)
            sp.write_bed(regs, peak_bed)
        n_peaks = sum(1 for _ in open(peak_bed))

        if not os.path.exists(se_tsv):
            st = se.call_superenhancers(peak_bed, bw, promoters, se_tsv,
                                        label=row["label"], binned=binned,
                                        binsize=sp.BIN)
        else:
            df = pd.read_csv(se_tsv, sep="\t")
            st = {
                "n_peaks_in": n_peaks,
                "n_peaks_after_promoter_filter": -1,
                "n_stitched": len(df),
                "n_super": int(df["is_super"].sum()),
                "cutoff_signal": float(df.loc[df.is_super == 1, "signal"].min())
                if df["is_super"].any() else float("nan"),
                "super_bp": int((df.loc[df.is_super == 1, "end"]
                                 - df.loc[df.is_super == 1, "start"]).sum()),
                "all_bp": int((df["end"] - df["start"]).sum()),
            }

        hg38_bed = os.path.join(SEDIR, cohort, f"{label}.hg38.bed")
        if not os.path.exists(hg38_bed):
            df = pd.read_csv(se_tsv, sep="\t")
            tmp = hg38_bed + ".native"
            # liftOver validates BED fields: 5 must be an integer and 7 an
            # integer, so the real signal travels inside the name field.
            df[["chrom", "start", "end"]].assign(
                name=[f"E{i}|{s:.4f}" for i, s in
                      enumerate(df["signal"].to_numpy(), 1)],
                score=0,
                strand=".",
                is_super=df["is_super"].astype(int),
            ).to_csv(tmp, sep="\t", header=False, index=False)
            if genome == "hg38":
                os.replace(tmp, hg38_bed)
            else:
                ok = liftover(tmp, hg38_bed,
                              os.path.join(REF, "hg19ToHg38.over.chain.gz"))
                os.remove(tmp)
                if not ok:
                    return {"label": row["label"], "cohort": cohort,
                            "status": "liftover-failed"}

        n_hg38 = sum(1 for _ in open(hg38_bed))
        return {
            "label": row["label"], "cohort": cohort, "gsm": row["gsm"],
            "genome": genome, "tissue_class": row["tissue_class"],
            "status": "ok", "n_peaks": n_peaks,
            "n_stitched": st["n_stitched"], "n_super": st["n_super"],
            "super_mb": st["super_bp"] / 1e6, "all_mb": st["all_bp"] / 1e6,
            "cutoff_signal": st["cutoff_signal"], "n_regions_hg38": n_hg38,
        }
    except Exception as exc:                      # keep the batch alive
        return {"label": row["label"], "cohort": cohort,
                "status": f"error: {type(exc).__name__}: {exc}"[:200]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--cohort", action="append")
    a = ap.parse_args()

    par, jac = params()
    print(f"peak-caller parameters: {par}" + (f"  (calibration J={jac:.3f})" if jac else ""))

    mf = pd.read_csv(os.path.join(DATA, "manifest.csv"))
    sel = mf[(mf["mark"] == "H3K27ac") & (mf["kind"] == "signal")].copy()
    if a.cohort:
        sel = sel[sel["cohort"].isin(a.cohort)]
    sel["signal_file"] = sel.apply(signal_path, axis=1)
    sel = sel[sel["signal_file"].map(os.path.exists)]
    print(f"{len(sel)} H3K27ac signal tracks ready")
    if sel.empty:
        return 0

    rows = []
    with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(one_sample, r, par) for _, r in sel.iterrows()]
        for i, f in enumerate(cf.as_completed(futs), 1):
            r = f.result()
            rows.append(r)
            msg = (f"SE {r.get('n_super')}/{r.get('n_stitched')} stitched"
                   if r["status"] == "ok" else r["status"])
            print(f"[{i}/{len(futs)}] {r['cohort']:<22} {str(r['label'])[:28]:<30} {msg}",
                  flush=True)

    out = pd.DataFrame(rows)
    path = os.path.join(RES, "se_summary.csv")
    if os.path.exists(path):
        old = pd.read_csv(path)
        key = ["cohort", "label"]
        out = pd.concat([old[~old.set_index(key).index.isin(out.set_index(key).index)],
                         out], ignore_index=True)
    out.to_csv(path, index=False)

    ok = out[out["status"] == "ok"]
    if not ok.empty:
        with pd.option_context("display.width", 200):
            print("\n== per-cohort super-enhancer counts ==")
            print(
                ok.groupby(["cohort", "tissue_class"])
                .agg(n=("label", "size"), se_med=("n_super", "median"),
                     se_iqr=("n_super", lambda s: s.quantile(.75) - s.quantile(.25)),
                     stitched_med=("n_stitched", "median"),
                     super_mb_med=("super_mb", "median"))
                .round(2).to_string()
            )
    bad = out[out["status"] != "ok"]
    if not bad.empty:
        print(f"\n{len(bad)} samples not done:")
        print(bad[["cohort", "label", "status"]].head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
