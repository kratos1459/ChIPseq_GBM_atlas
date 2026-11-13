#!/usr/bin/env bash
# Promoter windows on hg19, for the one cohort published on hg19. Peak calling
# and stitching happen in each cohort's native build; only the final enhancer
# coordinates are lifted to hg38 for cross-cohort comparison.
set -euo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chipseq
REF="$HOME/ChIPseq_GBM_atlas/data/ref"
cd "$REF"
# liftOver can place a TSS on an alt/random scaffold; keep main chromosomes only
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' gencode.v44.tss.hg19.bed \
  | sort -k1,1 -k2,2n > gencode.v44.tss.hg19.main.bed
bedtools slop -i gencode.v44.tss.hg19.main.bed -g hg19.main.chrom.sizes -b 2000 \
  | sort -k1,1 -k2,2n > promoters_2kb.hg19.bed
echo "hg19 promoters: $(wc -l < promoters_2kb.hg19.bed)"
