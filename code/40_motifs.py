#!/usr/bin/env python3
"""Transcription-factor motif enrichment in reproducible versus sporadic
super-enhancer domains.

The question is not "which motifs are in GBM enhancers" - that is answered by
any single dataset - but which motifs distinguish the domains that reproduce
across independent cohorts from those seen in a single sample. So the primary
set is always compared against a length-matched set of domains drawn from the
contrasting class, not against shuffled genomic background alone.

Comparisons
  core        vs sporadic           reproducible versus one-sample domains
  gsc_specific vs patient_specific  model-system versus tumour-specific
  core        vs shuffled           SEA's own control, as a sanity check

Outputs (results/motifs/)
  <comparison>/sea.tsv        raw SEA output
  motif_enrichment.csv        parsed, every comparison in one table
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REF = os.path.join(ROOT, "data", "ref")
CC = os.path.join(ROOT, "results", "cross_cohort")
OUT = os.path.join(ROOT, "results", "motifs")
FASTA = os.path.join(REF, "hg38.fa")
MOTIFS = os.path.join(REF, "jaspar_core_vertebrates.meme")
RNG = np.random.default_rng(20261113)
MAX_SEQ = 6000


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def length_match(primary, pool, rng=RNG):
    """Draw from pool so its length distribution matches primary's."""
    if not pool:
        return []
    plen = np.array([e - s for _, s, e in primary])
    pool_len = np.array([e - s for _, s, e in pool])
    edges = np.unique(np.quantile(plen, np.linspace(0, 1, 11)))
    want = np.histogram(plen, bins=edges)[0]
    chosen = []
    for i in range(len(want)):
        lo, hi = edges[i], edges[i + 1]
        cand = np.flatnonzero((pool_len >= lo) & (pool_len <= hi))
        if cand.size == 0:
            continue
        k = min(int(want[i]), cand.size)
        chosen += [pool[j] for j in rng.choice(cand, size=k, replace=False)]
    if len(chosen) < len(primary) // 2:      # fall back to a plain sample
        k = min(len(pool), len(primary))
        chosen = [pool[j] for j in rng.choice(len(pool), size=k, replace=False)]
    return chosen


def cap(regions, n=MAX_SEQ, rng=RNG):
    if len(regions) <= n:
        return regions
    return [regions[i] for i in rng.choice(len(regions), size=n, replace=False)]


def write_bed(regions, path):
    with open(path, "w") as fh:
        for i, (c, s, e) in enumerate(regions, 1):
            fh.write(f"{c}\t{s}\t{e}\tR{i}\n")


def getfasta(bed, out_fa):
    p = subprocess.run(
        ["bedtools", "getfasta", "-fi", FASTA, "-bed", bed, "-fo", out_fa],
        capture_output=True, text=True,
    )
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip()[:300])
    return out_fa


def run_sea(primary_fa, control_fa, outdir):
    os.makedirs(outdir, exist_ok=True)
    cmd = ["sea", "--p", primary_fa, "--m", MOTIFS, "--o", outdir, "--thresh", "1"]
    if control_fa:
        cmd[3:3] = ["--n", control_fa]
    p = subprocess.run(cmd, capture_output=True, text=True)
    tsv = os.path.join(outdir, "sea.tsv")
    if not os.path.exists(tsv):
        raise RuntimeError((p.stderr or p.stdout).strip()[-400:])
    return tsv


def parse_sea(tsv, comparison):
    rows = []
    with open(tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            d = dict(zip(header, f))
            rows.append(d)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["comparison"] = comparison
    for col in ("PVALUE", "QVALUE", "EVALUE", "ENR_RATIO", "TP%", "FP%"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def main():
    for need, what in ((FASTA, "hg38.fa"), (MOTIFS, "JASPAR motifs")):
        if not os.path.exists(need):
            print(f"missing {what} ({need}) - run 07_download_motif_ref.sh")
            return 1
    if shutil.which("sea") is None:
        print("MEME's sea not on PATH - activate the chipseq env")
        return 1
    os.makedirs(OUT, exist_ok=True)

    beds = {
        k: os.path.join(CC, f"{k}_domains.bed")
        for k in ("core", "sporadic", "gsc_specific", "patient_specific")
    }
    have = {k: v for k, v in beds.items() if os.path.exists(v)}
    print("domain sets:", {k: len(read_bed(v)) for k, v in have.items()})
    if "core" not in have:
        print("no core_domains.bed - run 30_cross_cohort.py first")
        return 1

    comparisons = [("core", "sporadic"), ("gsc_specific", "patient_specific"),
                   ("core", None)]
    frames = []
    for prim, ctrl in comparisons:
        if prim not in have or (ctrl and ctrl not in have):
            continue
        name = f"{prim}_vs_{ctrl or 'shuffled'}"
        wd = os.path.join(OUT, name)
        os.makedirs(wd, exist_ok=True)
        P = cap(read_bed(have[prim]))
        if len(P) < 20:
            print(f"{name}: only {len(P)} primary regions, skipped")
            continue
        pb, pf = os.path.join(wd, "primary.bed"), os.path.join(wd, "primary.fa")
        write_bed(P, pb)
        getfasta(pb, pf)
        cf_ = None
        if ctrl:
            C = length_match(P, read_bed(have[ctrl]))
            if len(C) < 20:
                print(f"{name}: control too small ({len(C)}), using shuffled")
            else:
                cb, cf_ = os.path.join(wd, "control.bed"), os.path.join(wd, "control.fa")
                write_bed(C, cb)
                getfasta(cb, cf_)
        print(f"{name}: {len(P)} primary"
              + (f", {len(C)} control" if cf_ else ", shuffled control"), flush=True)
        try:
            tsv = run_sea(pf, cf_, wd)
        except RuntimeError as exc:
            print(f"  SEA failed: {exc}")
            continue
        df = parse_sea(tsv, name)
        if df.empty:
            print("  SEA produced no rows")
            continue
        frames.append(df)
        col_q = "QVALUE" if "QVALUE" in df.columns else None
        sig = df[df[col_q] < 0.05] if col_q else df
        print(f"  {len(sig)} motifs at q<0.05 of {len(df)} tested")
        show = [c for c in ("ID", "ALT_ID", "ENR_RATIO", "QVALUE") if c in df.columns]
        if col_q:
            print(sig.nsmallest(12, col_q)[show].to_string(index=False))

    if frames:
        allm = pd.concat(frames, ignore_index=True)
        allm.to_csv(os.path.join(OUT, "motif_enrichment.csv"), index=False)
        print(f"\nwrote {len(allm)} rows -> results/motifs/motif_enrichment.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
