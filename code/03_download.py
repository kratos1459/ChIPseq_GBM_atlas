#!/usr/bin/env python3
"""Download the processed GEO files the atlas needs, plus reference annotation.

Groups (select with --group, repeatable, or --group all):
  small    GSE262089 HiChIP loops/peaks/RNA + GSE162976 narrowPeak/RNA  (~0.1 GB)
  probe    one GSE119776 bedGraph + one GSE145646 bigWig for format QC
  ref      hg19/hg38 chrom sizes, hg19->hg38 chain, GENCODE v44 basic, blacklists
  patient  GSE145646 H3K27ac bigWigs, 44 files                         (~9.0 GB)
  pedi     GSE162976 H3K27ac + H3K27M + H3K27me3 signal bigwigs        (~15.6 GB)
  gsc      GSE119776 H3K27ac bedGraphs, 98 files                       (~24.9 GB)

Everything resumes (wget -c) and is verified by size against GEO's manifest.
"""
import argparse
import pandas as pd
import concurrent.futures as cf
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RAW = os.path.join(DATA, "raw")
REF = os.path.join(DATA, "ref")


def ftp_sample_url(gse, gsm, fname):
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{fname}"
    )


REFS = [
    ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
     "hg38.chrom.sizes"),
    ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
     "hg19.chrom.sizes"),
    ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
     "hg19ToHg38.over.chain.gz"),
    ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
     "hg38ToHg19.over.chain.gz"),
    ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/"
     "gencode.v44.basic.annotation.gtf.gz", "gencode.v44.basic.annotation.gtf.gz"),
    ("https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/"
     "hg38-blacklist.v2.bed.gz", "hg38-blacklist.v2.bed.gz"),
    ("https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/"
     "hg19-blacklist.v2.bed.gz", "hg19-blacklist.v2.bed.gz"),
]


def wget(url, dest, expect=0):
    if os.path.exists(dest) and expect and os.path.getsize(dest) == expect:
        return dest, "cached", os.path.getsize(dest)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    cmd = ["wget", "-q", "-c", "--tries=8", "--timeout=120", "--waitretry=10",
           "-O", dest, url]
    p = subprocess.run(cmd, capture_output=True, text=True)
    got = os.path.getsize(dest) if os.path.exists(dest) else 0
    if p.returncode != 0:
        return dest, f"FAIL rc={p.returncode} {p.stderr.strip()[:120]}", got
    if expect and got != expect:
        return dest, f"SIZE-MISMATCH got={got} want={expect}", got
    return dest, "ok", got


def select(group, fl):
    """Return [(url, dest, expected_bytes)] for a group."""
    jobs = []
    if group == "small":
        sel = fl[
            ((fl.gse == "GSE262089") & fl.file.str.contains(r"\.bed(?:pe)?\.gz$|\.tsv\.gz$|\.txt\.gz$"))
            | ((fl.gse == "GSE162976") & fl.file.str.contains(r"narrowPeak\.gz$|tpm\.txt\.gz$"))
        ]
    elif group == "probe":
        a = fl[(fl.gse == "GSE119776") & fl.file.str.endswith("bedGraph.gz")].head(1)
        b = fl[(fl.gse == "GSE145646") & fl.file.str.endswith(".bw")].head(1)
        sel = pd.concat([a, b])
    elif group == "patient":
        sel = fl[(fl.gse == "GSE145646") & fl.file.str.endswith(".bw")]
    elif group == "pedi_k27ac":
        sel = fl[(fl.gse == "GSE162976")
                 & fl.file.str.endswith(".pval.signal.bigwig")
                 & fl.file.str.contains("H3K27AC", case=False)]
    elif group == "pedi":
        sel = fl[(fl.gse == "GSE162976") & fl.file.str.endswith(".pval.signal.bigwig")]
    elif group == "gsc":
        sel = fl[(fl.gse == "GSE119776") & fl.file.str.endswith(".bedGraph.gz")]
    else:
        raise SystemExit(f"unknown group {group}")
    for _, r in sel.iterrows():
        if not isinstance(r["gsm"], str):
            continue
        dest = os.path.join(RAW, r["gse"], r["file"])
        jobs.append((ftp_sample_url(r["gse"], r["gsm"], r["file"]), dest, int(r["bytes"])))
    return jobs


def main():
    import pandas as pd

    ap = argparse.ArgumentParser()
    ap.add_argument("--group", action="append", required=True)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    groups = args.group
    if "all" in groups:
        groups = ["ref", "small", "probe", "patient", "pedi", "gsc"]

    fl = pd.read_csv(os.path.join(DATA, "filelist.csv"))
    fl = fl[fl["gsm"].notna()]

    jobs = []
    for g in groups:
        if g == "ref":
            jobs += [(u, os.path.join(REF, n), 0) for u, n in REFS]
        else:
            jobs += select(g, fl)
    seen, uniq = set(), []          # groups can overlap; keep first occurrence
    for j in jobs:
        if j[1] in seen:
            continue
        seen.add(j[1])
        uniq.append(j)
    jobs = uniq
    want = sum(j[2] for j in jobs) / 1e9
    print(f"{len(jobs)} files, {want:.2f} GB expected -> {RAW}", flush=True)

    done = fail = 0
    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = {ex.submit(wget, *j): j for j in jobs}
        for i, fut in enumerate(cf.as_completed(futs), 1):
            dest, status, got = fut.result()
            if status in ("ok", "cached"):
                done += 1
            else:
                fail += 1
            if status != "cached" or i % 20 == 0:
                print(f"[{i}/{len(jobs)}] {status:<22} {got/1e6:9.1f} MB  "
                      f"{os.path.basename(dest)}", flush=True)
    print(f"\ndone={done} failed={fail}")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
