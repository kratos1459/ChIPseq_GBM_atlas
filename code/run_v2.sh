#!/usr/bin/env bash
# Re-run the downstream stages under the v2 definitions (patch_v2_definitions.py).
# Per-sample peak and super-enhancer calls are reused unchanged.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chipseq
ROOT="$HOME/ChIPseq_GBM_atlas"; LOG="$ROOT/logs"; cd "$ROOT" || exit 1
say() { echo "[$(date '+%F %T')] v2: $*" | tee -a "$LOG/pipeline.log"; }
say "start"
nice -n 5 python code/30_cross_cohort.py > "$LOG/v2_30_cross_cohort.log" 2>&1; say "30 exit $?"
nice -n 5 python code/40_motifs.py       > "$LOG/v2_40_motifs.log"       2>&1; say "40 exit $?"
nice -n 5 python code/50_expression.py   > "$LOG/v2_50_expression.log"   2>&1; say "50 exit $?"
nice -n 5 python code/60_figures.py      > "$LOG/v2_60_figures.log"      2>&1; say "60 exit $?"
/home/arian/miniconda3/bin/python3 code/70_write_report.py > "$LOG/v2_70_report.log" 2>&1; say "70 exit $?"
say "done"
