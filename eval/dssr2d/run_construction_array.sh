#!/bin/bash
# DSSR 2D of the construction-sweep predictions (needs a local DSSR binary; set DSSR=...). Usage:
#   mkdir -p results/analysis/logs/dssr2d; sbatch --array=0-9 eval/dssr2d/run_construction_array.sh
#   then: python eval/dssr2d/construction_2d.py --merge
#SBATCH -p cpu-single
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -J dssr2d
#SBATCH -o results/analysis/logs/dssr2d/dssr2d_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config.sh
NT=$(( ${SLURM_ARRAY_TASK_MAX:-0} + 1 ))
"$EVAL_PY" eval/dssr2d/construction_2d.py --task "${SLURM_ARRAY_TASK_ID}" --ntasks "$NT"
