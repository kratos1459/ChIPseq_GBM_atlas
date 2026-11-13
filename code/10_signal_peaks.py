#!/usr/bin/env python3
"""Call H3K27ac-enriched regions from a normalised signal track (bigWig).

Two of the four cohorts published signal tracks but no peak calls, so peaks have
to be derived from signal. Rather than pick a threshold by eye, this module is
built to be calibrated: 11_calibrate_caller.py tunes (quantile, gap, min_width)
against the one cohort that published BOTH signal and MACS2 peaks for the same
samples, and the tuned values are then applied to every cohort.

Method: bin the genome once per sample (cached to .npz), threshold bins at a
per-sample quantile of the non-zero signal distribution (so sequencing depth and
track scaling cancel), merge surviving bins across gaps, then drop short regions
and anything touching the ENCODE blacklist.
"""
import argparse
import gzip
import os

import numpy as np

BIN = 200


def read_sizes(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


def bin_genome(bigwig, chrom_sizes, binsize=BIN, cache_dir=None):
    """Mean signal per bin for every main chromosome. Returns dict chrom->array."""
    import pyBigWig

    cache = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        cache = os.path.join(
            cache_dir,
            os.path.basename(bigwig).replace(".", "_") + f".b{binsize}.npz",
        )
        if os.path.exists(cache):
            z = np.load(cache)
            return {k: z[k] for k in z.files}

    bw = pyBigWig.open(bigwig)
    try:
        have = bw.chroms()
        out = {}
        for c, n in chrom_sizes.items():
            if c not in have:
                continue
            length = min(n, have[c])
            nb = max(1, length // binsize)
            # Reading base-resolution values and averaging is ~25x faster than
            # bw.stats(nBins=...) and gives identical bin means; one chromosome
            # at a time keeps peak memory around 1 GB.
            a = bw.values(c, 0, nb * binsize, numpy=True)
            np.nan_to_num(a, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
            out[c] = a.reshape(nb, binsize).mean(axis=1).astype(np.float32)
            del a
    finally:
        bw.close()
    if cache:
        np.savez_compressed(cache, **out)
    return out


def nonzero_pool(binned, max_n=2_000_000):
    pool = []
    for a in binned.values():
        nz = a[a > 0]
        if nz.size:
            step = max(1, nz.size // max_n)
            pool.append(nz[::step])
    return np.concatenate(pool) if pool else np.array([], dtype=np.float32)


def call_from_binned(binned, quantile=0.99, gap=400, min_width=500, binsize=BIN,
                     pool=None, min_fold=None):
    """Return [(chrom, start, end, mean, max)] for one parameter set."""
    if pool is None:
        pool = nonzero_pool(binned)
    if pool.size == 0:
        return []
    thr = float(np.quantile(pool, quantile))
    if min_fold:
        thr = max(thr, min_fold * float(pool.mean()))

    gapbins = int(np.ceil(gap / binsize))
    regions = []
    for c, a in binned.items():
        idx = np.flatnonzero(a >= thr)
        if idx.size == 0:
            continue
        breaks = np.flatnonzero(np.diff(idx) > gapbins)
        starts = np.r_[idx[0], idx[breaks + 1]]
        ends = np.r_[idx[breaks], idx[-1]]
        for s, e in zip(starts, ends):
            bp_s, bp_e = int(s * binsize), int((e + 1) * binsize)
            if bp_e - bp_s < min_width:
                continue
            seg = a[s : e + 1]
            regions.append((c, bp_s, bp_e, float(seg.mean()), float(seg.max())))
    regions.sort(key=lambda r: (r[0], r[1]))
    return regions


def load_bed_intervals(path, main_only=True):
    """chrom -> sorted list of (start, end). Handles .gz."""
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if len(f) < 3:
                continue
            if main_only and ("_" in f[0] or f[0] in ("chrM", "chrEBV")):
                continue
            try:
                out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
            except ValueError:
                continue
    for c in out:
        out[c].sort()
    return out


def drop_blacklist(regions, blacklist_bed):
    bl = load_bed_intervals(blacklist_bed)
    keep = []
    for c, s, e, m, mx in regions:
        hit = False
        for bs, be in bl.get(c, ()):
            if bs >= e:
                break
            if s < be and bs < e:
                hit = True
                break
        if not hit:
            keep.append((c, s, e, m, mx))
    return keep


def write_bed(regions, path, name_prefix="R"):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w") as fh:
        for i, (c, s, e, m, mx) in enumerate(regions, 1):
            fh.write(f"{c}\t{s}\t{e}\t{name_prefix}_{i}\t{m:.5f}\t.\t{mx:.5f}\n")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bigwig")
    ap.add_argument("out_bed")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--quantile", type=float, default=0.99)
    ap.add_argument("--gap", type=int, default=400)
    ap.add_argument("--min-width", type=int, default=500)
    ap.add_argument("--blacklist")
    ap.add_argument("--cache-dir")
    a = ap.parse_args()

    sizes = read_sizes(a.chrom_sizes)
    binned = bin_genome(a.bigwig, sizes, cache_dir=a.cache_dir)
    regs = call_from_binned(binned, a.quantile, a.gap, a.min_width)
    if a.blacklist:
        regs = drop_blacklist(regs, a.blacklist)
    write_bed(regs, a.out_bed)
    tot = sum(e - s for _, s, e, _, _ in regs)
    print(f"{os.path.basename(a.bigwig)}: {len(regs)} regions, {tot/1e6:.1f} Mb")


if __name__ == "__main__":
    main()
