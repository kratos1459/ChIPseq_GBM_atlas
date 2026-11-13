#!/usr/bin/env bash
# Stages that run once every sample has been called: cross-cohort recurrence,
# motifs, expression, figures, report, local commit. Waits for run_pipeline.sh
# to finish stage C, then runs each stage in order. Resumable.
set -uo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chipseq

ROOT="$HOME/ChIPseq_GBM_atlas"
LOG="$ROOT/logs"
cd "$ROOT"
say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG/pipeline.log"; }

say "downstream: waiting for per-sample calling to finish"
for _ in $(seq 1 900); do
  if [ -s "$ROOT/results/se_summary.csv" ] \
     && ! pgrep -f "code/13_call_all.py" > /dev/null \
     && ! pgrep -f "code/run_pipeline.sh" > /dev/null; then
    break
  fi
  sleep 60
done
n_ok=$(awk -F, 'NR>1 && $0 ~ /,ok,/ {c++} END{print c+0}' "$ROOT/results/se_summary.csv" 2>/dev/null)
say "downstream: starting with ${n_ok:-0} completed samples"

run() {  # script log_name description
  say "downstream: $3"
  python "code/$1" > "$LOG/$2" 2>&1
  local rc=$?
  say "downstream: $3 exit $rc"
  tail -6 "$LOG/$2" | sed 's/^/    /' | tee -a "$LOG/pipeline.log"
}

run 30_cross_cohort.py  30_cross_cohort.log "cross-cohort recurrence and ceilings"
run 40_motifs.py        40_motifs.log       "motif enrichment"
run 50_expression.py    50_expression.log   "expression validation"
run 60_figures.py       60_figures.log      "figures"

# the report needs python-docx, which lives in the base interpreter
say "downstream: Methods and Results document"
python3 code/70_write_report.py > "$LOG/70_report.log" 2>&1
say "downstream: report exit $? -> $(tail -1 "$LOG/70_report.log")"

bash code/90_sync_github.sh "ChIP-seq GBM atlas: full cross-cohort results" \
  > "$LOG/90_sync.log" 2>&1
say "downstream: local commit -> $(tail -2 "$LOG/90_sync.log" | head -1)"
say "downstream: finished"
