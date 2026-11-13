#!/usr/bin/env bash
# GSE119776 published MACS2 pileup bedGraphs (hg19). Convert each to bigWig so
# every cohort is read through the same pyBigWig code path.
#
# The published files carry a UCSC track line, are not coordinate-sorted, and
# include scaffolds, so each is filtered to main chromosomes, sorted, clipped to
# the chromosome sizes and converted. Already-converted samples are skipped,
# so this can be re-run while downloads continue.
set -uo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chipseq

ROOT="$HOME/ChIPseq_GBM_atlas"
IN="$ROOT/data/raw/GSE119776"
OUT="$ROOT/data/bigwig_hg19"
SIZES="$ROOT/data/ref/hg19.main.chrom.sizes"
TMP="${TMPDIR:-/tmp}/bg2bw.$$"
JOBS="${JOBS:-3}"
mkdir -p "$OUT" "$TMP"
trap 'rm -rf "$TMP"' EXIT

convert() {
  local gz="$1" base out
  base="$(basename "$gz" .bedGraph.gz)"
  out="$OUT/$base.bw"
  [ -s "$out" ] && { echo "skip $base"; return 0; }
  local sorted="$TMP/$base.sorted.bg"
  # drop the track line and scaffolds, clamp to chromosome ends, sort
  zcat "$gz" \
    | awk -v OFS='\t' 'NR>1 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $2 ~ /^[0-9]+$/ {print $1,$2,$3,$4}' \
    | sort -k1,1 -k2,2n -S 1G -T "$TMP" > "$sorted"
  if [ ! -s "$sorted" ]; then echo "EMPTY $base"; rm -f "$sorted"; return 1; fi
  if bedGraphToBigWig "$sorted" "$SIZES" "$out.part" 2> "$TMP/$base.err"; then
    mv "$out.part" "$out"
    echo "ok $base $(du -h "$out" | cut -f1)"
  else
    echo "FAIL $base: $(tail -1 "$TMP/$base.err")"
    rm -f "$out.part"
  fi
  rm -f "$sorted"
}
export -f convert
export OUT SIZES TMP

shopt -s nullglob
files=("$IN"/*.bedGraph.gz)
echo "$(date) converting ${#files[@]} bedGraph files with $JOBS workers"
printf '%s\n' "${files[@]}" | xargs -P "$JOBS" -I{} bash -c 'convert "$@"' _ {}
echo "$(date) done: $(ls "$OUT"/*.bw 2>/dev/null | wc -l) bigWigs in $OUT"
