#!/usr/bin/env bash
# Resumable bulk download, calibration cohort first. Safe to re-run: files whose
# size already matches GEO's manifest are skipped.
ROOT="/home/arian/ChIPseq_GBM_atlas"
cd "$ROOT" || exit 1
exec nice -n 5 python3 "$ROOT/code/03_download.py" \
  --group pedi_k27ac --group patient --group pedi --group gsc --threads 4
