#!/usr/bin/env bash
# Independent expression data for validating the enhancer calls.
#
# TCGA GBM bulk RNA-seq from the UCSC Xena hub (gene symbols, log2(norm+1)),
# plus the matching clinical matrix. DepMap CRISPR dependency data is NOT
# fetched here: depmap.org serves a bot-verification page to non-browser
# clients, so that file has to be downloaded manually if it is wanted.
set -uo pipefail
REF="$HOME/ChIPseq_GBM_atlas/data/validation"
mkdir -p "$REF"
cd "$REF"
HUB="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download"

fetch() {
  [ -s "$2" ] && { echo "have $2"; return 0; }
  wget -q -c --tries=6 --timeout=180 -O "$2.part" "$1" && mv "$2.part" "$2" \
    || { echo "FAILED $1"; rm -f "$2.part"; return 1; }
  echo "got $2 ($(du -h "$2" | cut -f1))"
}

fetch "$HUB/TCGA.GBM.sampleMap%2FHiSeqV2.gz" tcga_gbm_hiseqv2.tsv.gz
fetch "$HUB/TCGA.GBM.sampleMap%2FGBM_clinicalMatrix" tcga_gbm_clinical.tsv
fetch "$HUB/TCGA.LGG.sampleMap%2FHiSeqV2.gz" tcga_lgg_hiseqv2.tsv.gz

for f in tcga_gbm_hiseqv2.tsv.gz tcga_lgg_hiseqv2.tsv.gz; do
  [ -s "$f" ] && echo "$f: $(zcat "$f" | head -1 | awk -F'\t' '{print NF-1}') samples, \
$(zcat "$f" | wc -l) rows"
done
echo "done $(date)"
