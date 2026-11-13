#!/usr/bin/env python3
"""One row per usable assay file: which GSM, cohort, tissue, mark, genome build,
and the local path to its signal track and/or peak calls.

Writes data/manifest.csv. Rows whose files are not on disk yet are kept with
present=False so the manifest is stable while downloads finish.
"""
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RAW = os.path.join(DATA, "raw")

GENOME = {
    "GSE119776": "hg19",
    "GSE145646": "hg38",
    "GSE162976": "hg38",
    "GSE262089": "hg38",
}
COHORT = {
    "GSE119776": "GSC_lines_adult",
    "GSE145646": "patient_tumour_adult",
    "GSE162976": "pHGG_pediatric",
    "GSE262089": "GSC_HiChIP",
}


def label_for(gse, title, fname):
    t = (title or "").strip()
    if gse == "GSE119776":
        return re.sub(r"_H3K27AC(_INPUT)?$", "", t, flags=re.I) or t
    if gse == "GSE145646":
        return re.sub(r"_H3K27ac$", "", t, flags=re.I) or t
    if gse == "GSE162976":
        m = re.match(r"(.+?)-(H3K27AC|H3K27M|H3K27ME3|ATAC|INPUT)", t, flags=re.I)
        base = m.group(1).strip() if m else t
        rep = re.search(r"Rep\s*(\d+)", t, flags=re.I)
        return f"{base}_rep{rep.group(1)}" if rep else base
    if gse == "GSE262089":
        m = re.match(r"(GSC|NSC)\s*(\w+)?-?(replicate(\d))?", t, flags=re.I)
        if m:
            kind = m.group(1).upper()
            ident = (m.group(2) or "").strip()
            rep = m.group(4)
            base = f"{kind}{ident}" if ident and ident.lower() != "replicate1" else kind
            return f"{base}_rep{rep}" if rep else base
    return t


SIGNAL_EXT = (".bw", ".bigwig", ".bedGraph.gz")
PEAK_PAT = re.compile(r"narrowPeak\.gz$|(?<!atac_)\w+\.bed\.gz$")


def kind_of(fname):
    f = fname.lower()
    if f.endswith(".bedpe.gz"):
        return "loops"
    if "narrowpeak" in f:
        return "peaks"
    if f.endswith(".bed.gz"):
        return "peaks_atac" if "atac_" in f else "peaks"
    if f.endswith((".bw", ".bigwig")):
        return "signal"
    if f.endswith(".bedgraph.gz"):
        return "signal"
    if f.endswith((".tpm.txt.gz", ".txt.gz", ".tsv.gz")):
        return "rna"
    return "other"


def main():
    sm = pd.read_csv(os.path.join(DATA, "sample_metadata.csv")).fillna("")
    fl = pd.read_csv(os.path.join(DATA, "filelist.csv"))
    fl = fl[fl["gsm"].notna()].copy()

    meta = sm.set_index(["gse", "gsm"]).to_dict("index")
    rows = []
    for _, r in fl.iterrows():
        key = (r["gse"], r["gsm"])
        m = meta.get(key, {})
        kind = kind_of(r["file"])
        if kind in ("other",):
            continue
        if r["file"].endswith(".idat.gz") or r["file"].endswith(".hic"):
            continue
        path = os.path.join(RAW, r["gse"], r["file"])
        mark = m.get("mark2") or m.get("mark") or ""
        if kind == "peaks_atac":
            mark = "ATAC-seq"
        rows.append(
            {
                "gse": r["gse"],
                "gsm": r["gsm"],
                "cohort": COHORT[r["gse"]],
                "genome": GENOME[r["gse"]],
                "label": label_for(r["gse"], m.get("title", ""), r["file"]),
                "title": m.get("title", ""),
                "mark": mark,
                "tissue_class": m.get("tissue_class", ""),
                "kind": "peaks" if kind == "peaks_atac" else kind,
                "file": r["file"],
                "path": path,
                "bytes": int(r["bytes"]),
                "present": os.path.exists(path)
                and os.path.getsize(path) == int(r["bytes"]),
            }
        )
    mf = pd.DataFrame(rows)

    # replicate-aware: is this an input/control track?
    mf["is_input"] = mf["mark"].eq("input")
    mf.to_csv(os.path.join(DATA, "manifest.csv"), index=False)

    with pd.option_context("display.width", 220, "display.max_rows", 120):
        print("== rows by cohort / kind / mark ==")
        print(
            mf.groupby(["cohort", "kind", "mark"])
            .agg(n=("gsm", "size"), present=("present", "sum"),
                 gb=("bytes", lambda s: round(s.sum() / 1e9, 2)))
            .to_string()
        )
        print("\n== H3K27ac signal tracks per cohort (the analysis backbone) ==")
        k = mf[(mf["mark"] == "H3K27ac") & (mf["kind"].isin(["signal", "peaks"]))]
        print(pd.crosstab([k["cohort"], k["kind"]], k["tissue_class"].replace("", "?")))
        print("\n== example labels ==")
        for c, d in k.groupby("cohort"):
            print(f"  {c}: {sorted(d['label'].unique())[:8]} ... n={d['label'].nunique()}")
    print(f"\nwrote {len(mf)} rows -> data/manifest.csv")


if __name__ == "__main__":
    main()
