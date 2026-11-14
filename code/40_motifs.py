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


PEAK_DIR = os.path.join(ROOT, "data", "peaks_called", "patient_tumour_adult")
HALF = 250


def peak_centres(domain_bed, workdir, tag):
    """Constituent H3K27ac peaks inside a domain set, as peak-centre +/-250 bp.

    SEA scores motif presence per sequence. Whole stitched domains are 10-40 kb,
    long enough to contain almost every motif, so primary and control both sit
    near 100% and nothing can be enriched. Testing the constituent peaks keeps
    the sequences short enough for presence to be informative. Peaks come from
    the hg38 patient-tumour calls, so no lifting is involved.
    """
    import glob as _g

    files = sorted(_g.glob(os.path.join(PEAK_DIR, "*.bed")))
    if not files:
        return []
    allp = os.path.join(workdir, f"_{tag}_allpeaks.bed")
    with open(allp, "w") as out:
        for f in files:
            with open(f) as fh:
                for line in fh:
                    p = line.split()
                    if len(p) >= 3:
                        out.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")
    merged = subprocess.run(
        f"sort -k1,1 -k2,2n {allp} | bedtools merge -i - | "
        f"bedtools intersect -u -a - -b {domain_bed}",
        shell=True, capture_output=True, text=True).stdout
    os.remove(allp)
    out = []
    for line in merged.splitlines():
        c, s, e = line.split()[:3]
        mid = (int(s) + int(e)) // 2
        out.append((c, max(0, mid - HALF), mid + HALF))
    return out


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
    # --text writes results to stdout. The bioconda MEME 5.5.9 build ships
    # without its HTML template, so the default mode writes only the TSV
    # header and then aborts ("Template does not contain data section").
    cmd = ["sea", "--text", "--p", primary_fa, "--m", MOTIFS, "--thresh", "1000"]
    if control_fa:
        cmd[2:2] = ["--n", control_fa]
    p = subprocess.run(cmd, capture_output=True, text=True)
    tsv = os.path.join(outdir, "sea.tsv")
    with open(tsv, "w") as fh:
        fh.write(p.stdout)
    if p.returncode != 0 or p.stdout.count("\n") < 2:
        raise RuntimeError((p.stderr or "no output").strip()[-400:])
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
        P = cap(peak_centres(have[prim], wd, "p"))
        if len(P) < 20:
            print(f"{name}: only {len(P)} primary regions, skipped")
            continue
        pb, pf = os.path.join(wd, "primary.bed"), os.path.join(wd, "primary.fa")
        write_bed(P, pb)
        getfasta(pb, pf)
        cf_ = None
        if ctrl:
            _pool = peak_centres(have[ctrl], wd, "c"); C = [_pool[i] for i in RNG.choice(len(_pool), size=min(len(_pool), len(P)), replace=False)] if _pool else []
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

    # base composition of every sequence set: a GC-box motif "enriched" against
    # a control of lower GC is composition, not transcription-factor binding
    import json as _json
    gc = {}
    for comp in os.listdir(OUT):
        for role in ("primary", "control"):
            fa = os.path.join(OUT, comp, f"{role}.fa")
            if os.path.exists(fa):
                seq = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
                if seq:
                    gc.setdefault(comp, {})[role] = round(
                        (seq.count("G") + seq.count("C")) / len(seq), 4)
    with open(os.path.join(OUT, "gc_composition.json"), "w") as fh:
        _json.dump(gc, fh, indent=2)
    print("GC composition:", gc)

    if frames:
        allm = pd.concat(frames, ignore_index=True)
        allm.to_csv(os.path.join(OUT, "motif_enrichment.csv"), index=False)
        print(f"\nwrote {len(allm)} rows -> results/motifs/motif_enrichment.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
