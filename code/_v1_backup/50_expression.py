#!/usr/bin/env python3
"""Does super-enhancer recurrence track transcription?

Each gene is given the maximum super-enhancer recurrence found within a window
around its TSS, separately for each cohort group. That per-gene recurrence is
then tested against expression measured in data that played no part in the
enhancer calls:

  * TCGA GBM bulk RNA-seq (n=172) and TCGA LGG (n=530) from the Xena hub;
  * pediatric pHGG RNA-seq TPM published alongside the ChIP (GSE162976);
  * GSC line RNA-seq counts from the HiChIP series (GSE262089).

Reported: Spearman correlation across all genes, and the expression of genes
carrying a reproducible ("core") super-enhancer against genes whose
super-enhancer was seen in a single sample, with Cliff's delta as effect size.

Outputs (results/expression/)
  gene_recurrence_expression.csv.gz   per-gene table
  expression_tests.csv                every statistical comparison
  summary.json
"""
import glob
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sstats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REF = os.path.join(DATA, "ref")
VAL = os.path.join(DATA, "validation")
RES = os.path.join(ROOT, "results")
CC = os.path.join(RES, "cross_cohort")
OUT = os.path.join(RES, "expression")
WINDOW = 100_000          # median HiChIP contact distance, from 20_loop_targets
BINSIZE = 5000


def cliffs_delta(a, b):
    """Rank-based effect size; +1 means a is entirely above b."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 5 or b.size < 5:
        return np.nan
    u = sstats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2 * u / (a.size * b.size) - 1)


def mw(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 5 or b.size < 5:
        return np.nan, np.nan
    r = sstats.mannwhitneyu(a, b, alternative="two-sided")
    return float(r.statistic), float(r.pvalue)


def load_genes():
    rows = []
    with open(os.path.join(REF, "gencode.v44.genes.hg38.bed")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, start, end, name, _, strand, gtype, gid = f[:8]
            tss = int(start) if strand == "+" else int(end) - 1
            rows.append((chrom, tss, name, gtype, gid.split(".")[0]))
    return pd.DataFrame(rows, columns=["chrom", "tss", "gene", "gene_type", "ensg"])


def per_gene_recurrence(genes, rec, window=WINDOW):
    """Max recurrence in any SE bin within +/- window of each TSS."""
    cols = [c for c in rec.columns if c.startswith("se_rec_")]
    out = {c: np.zeros(len(genes)) for c in cols}
    for chrom, gidx in genes.groupby("chrom").groups.items():
        sub = rec[rec["chrom"] == chrom]
        if sub.empty:
            continue
        starts = sub["start"].to_numpy()
        order = np.argsort(starts)
        starts = starts[order]
        vals = {c: sub[c].to_numpy()[order] for c in cols}
        gi = np.asarray(list(gidx))
        tss = genes.loc[gi, "tss"].to_numpy()
        lo = np.searchsorted(starts, tss - window, "left")
        hi = np.searchsorted(starts, tss + window, "right")
        for c in cols:
            v = vals[c]
            for a, b, pos in zip(lo, hi, gi):   # genes is reset_index'd: index == position
                out[c][pos] = v[a:b].max() if b > a else 0.0
    for c in cols:
        genes[c.replace("se_rec_", "rec_")] = out[c]
    return genes


def load_tcga(path):
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    return df.median(axis=1).rename("expr")


def load_geo_rna():
    """Pediatric TPM and GSC counts, averaged per gene."""
    out = {}
    ped = sorted(glob.glob(os.path.join(DATA, "raw", "GSE162976", "*tpm.txt.gz")))
    if ped:
        frames = []
        for p in ped:
            d = pd.read_csv(p, sep="\t")
            if d.shape[1] < 2:
                continue
            d.columns = ["gene_id", "tpm"][: d.shape[1]]
            d["ensg"] = d["gene_id"].astype(str).str.split(".").str[0]
            frames.append(d.groupby("ensg")["tpm"].max())
        if frames:
            m = pd.concat(frames, axis=1)
            out["pedi_tpm"] = np.log2(m.median(axis=1) + 1).rename("expr")
            print(f"  pediatric RNA: {len(frames)} samples, {len(out['pedi_tpm'])} genes")
    gsc = [p for p in sorted(glob.glob(os.path.join(DATA, "raw", "GSE262089", "*.txt.gz")))
           if "barcodes" not in p and "features" not in p]
    if gsc:
        frames = []
        for p in gsc:
            d = pd.read_csv(p, sep="\t", header=None, names=["ensg", "count"])
            d["ensg"] = d["ensg"].astype(str).str.split(".").str[0]
            s = d.groupby("ensg")["count"].sum()
            frames.append(1e6 * s / s.sum())            # CPM
        m = pd.concat(frames, axis=1)
        out["gsc_cpm"] = np.log2(m.median(axis=1) + 1).rename("expr")
        print(f"  GSC RNA: {len(frames)} samples, {len(out['gsc_cpm'])} genes")
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    rpath = os.path.join(CC, "recurrence.csv.gz")
    if not os.path.exists(rpath):
        print("no recurrence.csv.gz - run 30_cross_cohort.py first")
        return 1
    rec = pd.read_csv(rpath)
    print(f"{len(rec)} super-enhancer bins; groups: "
          f"{[c[len('se_rec_'):] for c in rec.columns if c.startswith('se_rec_')]}")

    genes = load_genes().reset_index(drop=True)
    genes = per_gene_recurrence(genes, rec)
    print(f"{len(genes)} genes annotated with recurrence in +/-{WINDOW//1000} kb")

    expr = {}
    for name, f in (("tcga_gbm", "tcga_gbm_hiseqv2.tsv.gz"),
                    ("tcga_lgg", "tcga_lgg_hiseqv2.tsv.gz")):
        p = os.path.join(VAL, f)
        if os.path.exists(p):
            expr[name] = load_tcga(p)
            print(f"  {name}: {len(expr[name])} genes")
    expr.update(load_geo_rna())

    g = genes.copy()
    for name, s in expr.items():
        if name.startswith("tcga"):
            g[name] = g["gene"].map(s)
        else:
            g[name] = g["ensg"].map(s)

    has_gsc = "rec_GSC_line" in g.columns
    has_pat = "rec_patient_tumour" in g.columns
    if has_gsc and has_pat:
        g["rec_adult_min"] = g[["rec_GSC_line", "rec_patient_tumour"]].min(axis=1)
        g["rec_adult_mean"] = g[["rec_GSC_line", "rec_patient_tumour"]].mean(axis=1)
        nsc = g["rec_NSC"] if "rec_NSC" in g.columns else 0.0
        g["class"] = np.where(
            (g["rec_adult_min"] >= 0.5) & (nsc < 0.10), "core",
            np.where(g["rec_adult_mean"] > 0, "other_se", "no_se"),
        )
    g.to_csv(os.path.join(OUT, "gene_recurrence_expression.csv.gz"),
             index=False, compression="gzip")

    rows = []
    pc = [c for c in g.columns if c.startswith("rec_")]
    for ename in expr:
        sub = g[np.isfinite(g[ename])]
        for rcol in pc:
            if sub[rcol].nunique() < 3:
                continue
            r = sstats.spearmanr(sub[rcol], sub[ename])
            rows.append({"test": "spearman", "expression": ename,
                         "recurrence": rcol, "n": int(len(sub)),
                         "rho": float(r.statistic), "p": float(r.pvalue)})
        if "class" in sub.columns:
            core = sub.loc[sub["class"] == "core", ename]
            other = sub.loc[sub["class"] == "other_se", ename]
            none = sub.loc[sub["class"] == "no_se", ename]
            for a, b, lab in ((core, other, "core_vs_other_se"),
                              (core, none, "core_vs_no_se"),
                              (other, none, "other_se_vs_no_se")):
                u, p = mw(a, b)
                rows.append({
                    "test": lab, "expression": ename, "recurrence": "class",
                    "n": int(len(a)), "n_b": int(len(b)),
                    "median_a": float(np.nanmedian(a)) if len(a) else np.nan,
                    "median_b": float(np.nanmedian(b)) if len(b) else np.nan,
                    "U": u, "p": p, "cliffs_delta": cliffs_delta(a, b),
                })
    if "tcga_gbm" in expr and "tcga_lgg" in expr and "class" in g.columns:
        g["gbm_minus_lgg"] = g["tcga_gbm"] - g["tcga_lgg"]
        core = g.loc[g["class"] == "core", "gbm_minus_lgg"]
        other = g.loc[g["class"] == "other_se", "gbm_minus_lgg"]
        u, p = mw(core, other)
        rows.append({"test": "gbm_vs_lgg_delta_core_vs_other", "expression":
                     "tcga_gbm-lgg", "recurrence": "class", "n": int(core.notna().sum()),
                     "n_b": int(other.notna().sum()),
                     "median_a": float(np.nanmedian(core)),
                     "median_b": float(np.nanmedian(other)), "U": u, "p": p,
                     "cliffs_delta": cliffs_delta(core, other)})

    tests = pd.DataFrame(rows)
    tests.to_csv(os.path.join(OUT, "expression_tests.csv"), index=False)
    with pd.option_context("display.width", 240, "display.max_rows", 80):
        print("\n== recurrence vs expression ==")
        print(tests.round(4).to_string(index=False))

    summary = {"window_kb": WINDOW // 1000, "n_genes": int(len(g))}
    if "class" in g.columns:
        summary["class_counts"] = g["class"].value_counts().to_dict()
        named = g.loc[(g["class"] == "core") & (~g["gene"].str.startswith("ENSG"))]
        if "tcga_gbm" in g.columns:
            top = named.nlargest(30, "tcga_gbm")["gene"].tolist()
            summary["top_core_genes_by_tcga_expression"] = top
            print("\ncore-SE genes most highly expressed in TCGA GBM:", top)
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print("\nsummary:", json.dumps(summary, indent=2, default=str)[:1500])
    return 0


if __name__ == "__main__":
    sys.exit(main())
