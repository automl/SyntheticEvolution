#!/bin/bash
# Regenerate every LaTeX table, assemble the copy-paste folder for the manuscript, and verify all numbers independently.
#   sbatch --export=ALL,MANUSCRIPT=<paper.tex>,DEST=<dir> tables/paper_tables.sh     (or run with bash)
# Needs results/dssr2d/construction_2d.tsv.gz (eval/dssr2d) for the construction-stage table.
#SBATCH -p cpu-single
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J paper_tables
#SBATCH -o results/analysis/logs/paper_tables_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
source config.sh
: "${MANUSCRIPT:?set MANUSCRIPT}" "${DEST:?set DEST}"
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
st() { echo "$1" | tee -a "$DEST/STATUS.md"; }
mkdir -p "$DEST"; echo "# Status of the CRC table build ($(date '+%F %T'), job ${SLURM_JOB_ID:-local})" > "$DEST/STATUS.md"
"$EVAL_PY" tables/make_crc_tables.py --min-preds 20 > /dev/null && st "- ok   make_crc_tables" || { st "- FAIL make_crc_tables"; exit 1; }
"$EVAL_PY" tables/make_extra_tables.py > /dev/null && st "- ok   make_extra_tables" || { st "- FAIL make_extra_tables"; exit 1; }
"$EVAL_PY" tables/prose_numbers.py > /dev/null && st "- ok   prose_numbers (results/tables/prose_numbers.md)" || { st "- FAIL prose_numbers"; exit 1; }
"$EVAL_PY" tables/assemble_paper_tables.py --manuscript "$MANUSCRIPT" --dest "$DEST" > /dev/null && st "- ok   assemble_paper_tables -> $DEST" || { st "- FAIL assemble"; exit 1; }
if "$EVAL_PY" tables/verify_tables.py --dest "$DEST" > /dev/null; then st "- ok   verify_tables: $(grep -o '[0-9]* checks, [0-9]* failed' "$DEST/VERIFY_REPORT.md")"
else st "- FAIL verify_tables: $(grep -o '[0-9]* checks, [0-9]* failed' "$DEST/VERIFY_REPORT.md") (see VERIFY_REPORT.md)"; exit 1; fi
st "- DONE"
