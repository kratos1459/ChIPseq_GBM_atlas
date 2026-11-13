#!/usr/bin/env bash
# Reference annotation for the atlas: main-chromosome sizes, TSS and gene body
# BEDs from GENCODE v44 (hg38), blacklists, and an hg19 projection of the TSS
# set so the hg19 cohort can be queried without lifting whole signal tracks.
set -euo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chipseq

ROOT="$HOME/ChIPseq_GBM_atlas"
REF="$ROOT/data/ref"
cd "$REF"

MAIN='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'

# --- main-chromosome size files -------------------------------------------
for g in hg38 hg19; do
  awk -v re="$MAIN" '$1 ~ re {print $1"\t"$2}' "$g.chrom.sizes" \
    | sort -k1,1 > "$g.main.chrom.sizes"
done
echo "chrom sizes: hg38=$(wc -l < hg38.main.chrom.sizes) hg19=$(wc -l < hg19.main.chrom.sizes)"

# --- blacklists -----------------------------------------------------------
for g in hg38 hg19; do
  zcat "$g-blacklist.v2.bed.gz" | cut -f1-3 | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "$g.blacklist.bed"
done
echo "blacklist: hg38=$(wc -l < hg38.blacklist.bed) regions"

# --- genes and TSS from GENCODE v44 ---------------------------------------
# keep protein_coding + lncRNA, level-verified gene records on main chroms
zcat gencode.v44.basic.annotation.gtf.gz \
 | awk -F'\t' -v re="$MAIN" '$3=="gene" && $1 ~ re {print}' \
 | python3 "$ROOT/code/04b_parse_gencode.py" \
       gencode.v44.genes.hg38.bed gencode.v44.tss.hg38.bed
for b in gencode.v44.genes.hg38.bed gencode.v44.tss.hg38.bed; do
  sort -k1,1 -k2,2n "$b" -o "$b"
done
echo "genes=$(wc -l < gencode.v44.genes.hg38.bed) tss=$(wc -l < gencode.v44.tss.hg38.bed)"
awk '{print $7}' gencode.v44.genes.hg38.bed | sort | uniq -c

# --- promoter windows (TSS +/- 2 kb), used to exclude promoters from SEs ---
bedtools slop -i gencode.v44.tss.hg38.bed -g hg38.main.chrom.sizes -b 2000 \
  | sort -k1,1 -k2,2n > promoters_2kb.hg38.bed
echo "promoters=$(wc -l < promoters_2kb.hg38.bed)"

# --- hg19 projection of the TSS set (for the hg19 GSC cohort) --------------
# liftOver demands an integer in field 7, so lift a BED4 keyed on gene_id and
# rejoin the annotation afterwards.
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8}' gencode.v44.tss.hg38.bed > _tss4.hg38.bed
liftOver _tss4.hg38.bed hg38ToHg19.over.chain.gz \
  _tss4.hg19.bed tss.hg38tohg19.unmapped 2> liftover_tss.log || true
awk 'BEGIN{OFS="\t"} NR==FNR{nm[$8]=$4; st[$8]=$6; ty[$8]=$7; next}
     ($4 in nm){print $1,$2,$3,nm[$4],0,st[$4],ty[$4],$4}' \
  gencode.v44.tss.hg38.bed _tss4.hg19.bed \
  | sort -k1,1 -k2,2n > gencode.v44.tss.hg19.bed
rm -f _tss4.hg38.bed _tss4.hg19.bed
echo "tss lifted to hg19: $(wc -l < gencode.v44.tss.hg19.bed) (unmapped $(grep -vc '^#' tss.hg38tohg19.unmapped || echo 0))"
echo "=== annotation done $(date) ==="
