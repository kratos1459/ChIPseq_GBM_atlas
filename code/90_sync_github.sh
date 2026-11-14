#!/usr/bin/env bash
# Mirror code and small result tables into ~/GitHub_ready/ChIPseq_GBM_atlas and
# keep a local git history there. Raw signal tracks, bigWigs, genome FASTA and
# binned caches are deliberately excluded - only things worth reviewing go in.
set -uo pipefail
SRC="$HOME/ChIPseq_GBM_atlas"
DST="$HOME/GitHub_ready/ChIPseq_GBM_atlas"
mkdir -p "$DST/code" "$DST/data" "$DST/results" "$DST/figures"

rsync -a --delete "$SRC/code/" "$DST/code/"
cp -f "$SRC/README.md" "$DST/README.md" 2>/dev/null
cp -f "$SRC/00_dataset_inventory.md" "$DST/" 2>/dev/null

# manifests and reference tables that are small and worth versioning
for f in sample_metadata.csv series_metadata.csv filelist.csv \
         filelist_annotated.csv suppl_inventory.csv manifest.csv; do
  [ -f "$SRC/data/$f" ] && cp -f "$SRC/data/$f" "$DST/data/$f"
done

# result tables and beds, excluding per-sample bulk
rsync -a --delete --delete-excluded \
  --exclude='se/' --exclude='se_test/' --exclude='_v1_*/' \
  --include='*/' \
  --include='*.csv' --include='*.csv.gz' --include='*.json' \
  --include='*.bed' --include='*.tsv' \
  --exclude='se/**' --exclude='*.fa' --exclude='*.fasta' \
  --exclude='*' \
  "$SRC/results/" "$DST/results/" 2>/dev/null
rsync -a --delete "$SRC/figures/" "$DST/figures/" 2>/dev/null
mkdir -p "$DST/manuscript"
rsync -a --delete "$SRC/manuscript/" "$DST/manuscript/" 2>/dev/null

cat > "$DST/.gitignore" <<'EOF'
# raw and derived sequencing data are too large to version
data/raw/
data/bigwig_hg19/
data/cache_bins/
data/ref/*.fa
data/ref/*.fa.gz
data/ref/*.chain.gz
data/ref/*.gtf.gz
data/validation/
results/se/
*.bw
*.bigwig
*.bedGraph.gz
*.npz
__pycache__/
*.pyc
EOF

cd "$DST" || exit 1
if [ ! -d .git ]; then
  git init -q
  git config user.name  "${GIT_NAME:-Tawsif Al Arian}"
  git config user.email "${GIT_EMAIL:-tawsifalarian43@gmail.com}"
fi
git add -A
if git diff --cached --quiet; then
  echo "nothing new to commit"
else
  git commit -q -m "${1:-ChIP-seq GBM enhancer atlas: pipeline and results update}

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01UNW4MtE7sGywWyw3jNmy8x"
  echo "committed: $(git log --oneline -1)"
fi
echo "files tracked: $(git ls-files | wc -l)  size: $(du -sh "$DST" | cut -f1)"
if git remote | grep -q origin; then
  git push -q origin HEAD && echo "pushed to origin"
else
  echo "no 'origin' remote yet - create the GitHub repository, then:"
  echo "  git -C $DST remote add origin <url> && git -C $DST push -u origin master"
fi
