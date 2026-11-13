#!/usr/bin/env bash
# Build the ChIP-seq analysis environment. Two stages so a solver failure in the
# optional stage cannot block the core tools.
set -u
source "$HOME/miniconda3/etc/profile.d/conda.sh"
export CONDA_SOLVER=libmamba

CORE="python=3.11 samtools bedtools macs3 deeptools bedops pybigwig pybedtools
      requests pandas scipy numpy statsmodels scikit-learn matplotlib seaborn
      pyarrow tqdm"
EXTRA="ucsc-bigwigaverageoverbed ucsc-liftover ucsc-bedgraphtobigwig
       ucsc-bigwigmerge meme sra-tools bowtie2 bedtools"

echo "=== $(date) stage 1: core ==="
conda create -y -n chipseq -c conda-forge -c bioconda $CORE 2>&1 | tail -25
echo "=== $(date) stage 1 exit: $? ==="

echo "=== $(date) stage 2: extras ==="
conda install -y -n chipseq -c conda-forge -c bioconda $EXTRA 2>&1 | tail -25
echo "=== $(date) stage 2 exit: $? ==="

conda activate chipseq
echo "=== verification ==="
for t in samtools bedtools macs3 bamCoverage computeMatrix plotHeatmap \
         bedmap bigWigAverageOverBed liftOver bedGraphToBigWig sea ame \
         fasterq-dump bowtie2; do
  printf "%-26s %s\n" "$t" "$(command -v "$t" || echo MISSING)"
done
python - <<'PY'
import importlib
for m in ["pyBigWig","pybedtools","pandas","scipy","sklearn","statsmodels","matplotlib","seaborn"]:
    try:
        importlib.import_module(m); print(m, "OK")
    except Exception as e:
        print(m, "MISSING", e)
PY
echo "=== $(date) done ==="
