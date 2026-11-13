#!/usr/bin/env bash
# Unattended pipeline: calibrate the peak caller as soon as the calibration
# cohort has landed, convert the hg19 bedGraph cohort once its downloads finish,
# then call enhancers and super-enhancers for every sample.
#
# Each stage is resumable and logs separately, so this can be killed and
# restarted without losing completed work.
set -uo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chipseq

ROOT="$HOME/ChIPseq_GBM_atlas"
LOG="$ROOT/logs"
cd "$ROOT"
say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG/pipeline.log"; }

# ---------------------------------------------------------------- stage A ----
# 21 pediatric H3K27ac bigwigs are the calibration set.
say "stage A: waiting for the calibration cohort (21 bigwigs)"
for _ in $(seq 1 240); do
  n=$(ls "$ROOT"/data/raw/GSE162976/*H3K27AC*.pval.signal.bigwig 2>/dev/null | wc -l)
  [ "$n" -ge 21 ] && break
  sleep 60
done
say "stage A: $n/21 calibration bigwigs present"
python3 code/05_build_manifest.py > "$LOG/05_manifest.log" 2>&1
python  code/11_calibrate_caller.py > "$LOG/11_calibrate.log" 2>&1
say "stage A: calibration exit $? -> $(grep -c . "$LOG/11_calibrate.log") log lines"
grep -A12 '^selected:' "$LOG/11_calibrate.log" | tee -a "$LOG/pipeline.log"

# ---------------------------------------------------------------- stage B ----
say "stage B: waiting for downloads to finish"
for _ in $(seq 1 480); do
  pgrep -f "code/03_download.py" > /dev/null || break
  sleep 60
done
python3 code/03b_verify_downloads.py > "$LOG/03b_verify.log" 2>&1
say "stage B: verification -> $(head -1 "$LOG/03b_verify.log")"

say "stage B: converting hg19 bedGraphs to bigWig"
JOBS=4 bash code/06_bedgraph_to_bigwig.sh > "$LOG/06_bedgraph2bw.log" 2>&1
say "stage B: $(ls "$ROOT"/data/bigwig_hg19/*.bw 2>/dev/null | wc -l) bigWigs built"

# ---------------------------------------------------------------- stage C ----
say "stage C: calling enhancers and super-enhancers for every sample"
python3 code/05_build_manifest.py > "$LOG/05_manifest.log" 2>&1
python  code/13_call_all.py --workers 4 > "$LOG/13_call_all.log" 2>&1
say "stage C: exit $?"
tail -24 "$LOG/13_call_all.log" | tee -a "$LOG/pipeline.log"
say "pipeline finished"
