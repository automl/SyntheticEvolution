#!/bin/bash
# SPOT-RNA topologies as a CPU array (see predict_2d.py). Usage:
#   mkdir -p results/analysis/logs/topologies; sbatch --array=0-19 analysis/topologies/run_spotrna_array.sh single   # then --merge
#   sbatch --array=0-19 analysis/topologies/run_spotrna_array.sh complex1v1
#SBATCH -p cpu-single
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 01:00:00
#SBATCH -J spotrna2d
#SBATCH -o results/analysis/logs/topologies/spotrna2d_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config.sh
NT=$(( ${SLURM_ARRAY_TASK_MAX:-0} + 1 ))
export TMPDIR=$(mktemp -d "${CRC_ROOT}/data/topologies/shards/tmp.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT
"$HOME/.conda/envs/${RNABENCH_ENV}/bin/python" analysis/topologies/predict_2d.py --method spotrna --set "$1" \
    --task "${SLURM_ARRAY_TASK_ID}" --ntasks "$NT"
