#!/usr/bin/env python3
"""ROSE-style stitching and super-enhancer ranking, reimplemented.

Steps, following Whyte/Loven 2013:
  1. drop peaks lying entirely within a promoter window (TSS +/- 2 kb), so that
     highly acetylated promoters are not mistaken for enhancer domains;
  2. stitch surviving peaks within 12.5 kb of each other;
  3. score every stitched region by the total H3K27ac signal it contains;
  4. rank regions by score and split them at the inflection point of the ranked
     curve - the point where the tangent to the normalised curve has slope 1.

The inflection is found numerically rather than by a fixed rank or percentile,
so the super-enhancer count is a property of each sample's own distribution.
"""
import argparse
import os

import numpy as np

STITCH = 12500
TSS_PAD = 2000


def read_peaks(path):
    import gzip

    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if len(f) < 3 or "_" in f[0] or f[0] in ("chrM", "chrEBV"):
                continue
            try:
                rows.append((f[0], int(f[1]), int(f[2])))
            except ValueError:
                continue
    rows.sort()
    return rows


def load_promoters(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                d.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    for c in d:
        d[c].sort()
    return d


def inside_any(ivs, s, e):
    """True if [s,e) lies entirely inside one of the sorted intervals."""
    if not ivs:
        return False
    import bisect

    starts = [a for a, _ in ivs]
    i = bisect.bisect_right(starts, s) - 1
    for j in (i, i + 1):
        if 0 <= j < len(ivs):
            a, b = ivs[j]
            if a <= s and e <= b:
                return True
    return False


def stitch(peaks, distance=STITCH):
    out = []
    cur_c = cur_s = cur_e = None
    n = 0
    for c, s, e in peaks:
        if cur_c == c and s - cur_e <= distance:
            cur_e = max(cur_e, e)
            n += 1
        else:
            if cur_c is not None:
                out.append((cur_c, cur_s, cur_e, n))
            cur_c, cur_s, cur_e, n = c, s, e, 1
    if cur_c is not None:
        out.append((cur_c, cur_s, cur_e, n))
    return out


def score_regions_binned(regions, binned, binsize=200):
    """Total signal per region, read off an already-binned genome."""
    scores = []
    for c, s, e, _ in regions:
        a = binned.get(c)
        if a is None:
            scores.append(0.0)
            continue
        i, j = s // binsize, min(len(a), (e + binsize - 1) // binsize)
        scores.append(float(a[i:j].sum() * binsize) if j > i else 0.0)
    return np.array(scores, dtype=float)


def score_regions(regions, bigwig):
    """Total signal (area) per stitched region."""
    import pyBigWig

    bw = pyBigWig.open(bigwig)
    try:
        have = bw.chroms()
        scores = []
        for c, s, e, _ in regions:
            if c not in have:
                scores.append(0.0)
                continue
            e2 = min(e, have[c])
            if e2 <= s:
                scores.append(0.0)
                continue
            try:
                v = bw.stats(c, s, e2, type="mean", exact=True)[0]
            except RuntimeError:
                v = None
            m = 0.0 if v is None or not np.isfinite(v) else float(v)
            scores.append(m * (e2 - s))
    finally:
        bw.close()
    return np.array(scores, dtype=float)


def inflection_cutoff(scores):
    """Index-free cutoff: value where the normalised ranked curve has slope 1."""
    y = np.sort(scores)                      # ascending
    if y.size < 10 or y.max() <= 0:
        return float("inf"), 0
    x = np.arange(y.size, dtype=float) / (y.size - 1)
    yn = (y - y.min()) / (y.max() - y.min())
    # smooth derivative over a window scaled to the data size
    w = max(5, y.size // 200)
    dy = np.gradient(yn, x, edge_order=2)
    k = np.ones(w) / w
    dys = np.convolve(dy, k, mode="same")
    below = np.flatnonzero(dys < 1.0)
    idx = int(below[-1]) if below.size else 0
    cutoff = float(y[idx])
    n_super = int((scores > cutoff).sum())
    return cutoff, n_super


def call_superenhancers(peaks_bed, bigwig, promoters_bed, out_tsv,
                        stitch_distance=STITCH, label=None, binned=None,
                        binsize=200):
    peaks = read_peaks(peaks_bed)
    prom = load_promoters(promoters_bed)
    kept = [(c, s, e) for c, s, e in peaks if not inside_any(prom.get(c, []), s, e)]
    regions = stitch(kept, stitch_distance)
    scores = (score_regions_binned(regions, binned, binsize) if binned is not None
              else score_regions(regions, bigwig))
    cutoff, n_super = inflection_cutoff(scores)

    order = np.argsort(-scores)
    d = os.path.dirname(out_tsv)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(out_tsv, "w") as fh:
        fh.write("chrom\tstart\tend\tn_peaks\tsignal\trank\tis_super\n")
        for rank, i in enumerate(order, 1):
            c, s, e, npk = regions[i]
            fh.write(
                f"{c}\t{s}\t{e}\t{npk}\t{scores[i]:.4f}\t{rank}\t"
                f"{int(scores[i] > cutoff)}\n"
            )
    stats = {
        "label": label or os.path.basename(peaks_bed),
        "n_peaks_in": len(peaks),
        "n_peaks_after_promoter_filter": len(kept),
        "n_stitched": len(regions),
        "n_super": n_super,
        "cutoff_signal": cutoff,
        "super_bp": int(sum(regions[i][2] - regions[i][1]
                            for i in range(len(regions)) if scores[i] > cutoff)),
        "all_bp": int(sum(e - s for _, s, e, _ in regions)),
    }
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("peaks")
    ap.add_argument("bigwig")
    ap.add_argument("out_tsv")
    ap.add_argument("--promoters", required=True)
    ap.add_argument("--stitch", type=int, default=STITCH)
    ap.add_argument("--label")
    a = ap.parse_args()
    st = call_superenhancers(a.peaks, a.bigwig, a.promoters, a.out_tsv,
                             a.stitch, a.label)
    print(
        f"{st['label']}: peaks {st['n_peaks_in']} -> {st['n_peaks_after_promoter_filter']} "
        f"-> stitched {st['n_stitched']}, SE {st['n_super']} "
        f"({st['super_bp']/1e6:.1f} Mb of {st['all_bp']/1e6:.1f} Mb)"
    )


if __name__ == "__main__":
    main()
