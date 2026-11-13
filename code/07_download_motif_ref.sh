#!/usr/bin/env bash
# Reference data for motif enrichment: hg38 sequence and the JASPAR core
# vertebrate motif set in MEME format. Both are large-ish and only needed for
# the motif stage, so they download separately from the signal tracks.
set -uo pipefail
REF="$HOME/ChIPseq_GBM_atlas/data/ref"
mkdir -p "$REF"
cd "$REF"

fetch() {  # url dest
  [ -s "$2" ] && { echo "have $2"; return 0; }
  echo "fetching $(basename "$2")"
  wget -q -c --tries=6 --timeout=180 -O "$2.part" "$1" && mv "$2.part" "$2" \
    || { echo "FAILED $1"; rm -f "$2.part"; return 1; }
}

fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" hg38.fa.gz

# JASPAR mirrors have moved more than once; try each in turn.
for u in \
  "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt" \
  "https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt" \
  "https://jaspar2020.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt"
do
  fetch "$u" jaspar_core_vertebrates.meme && break
done

if [ -s hg38.fa.gz ] && [ ! -s hg38.fa ]; then
  echo "decompressing hg38.fa.gz"
  gunzip -c hg38.fa.gz > hg38.fa.part && mv hg38.fa.part hg38.fa
fi
ls -la hg38.fa hg38.fa.gz jaspar_core_vertebrates.meme 2>/dev/null
grep -c '^MOTIF' jaspar_core_vertebrates.meme 2>/dev/null | sed 's/^/motifs: /'
echo "done $(date)"
