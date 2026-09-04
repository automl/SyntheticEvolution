#!/bin/bash
#SBATCH --job-name=trrosetta_shs
#SBATCH --partition=gpu-single
#SBATCH --output=logs/trrosetta_shs_%j.out
#SBATCH --error=logs/trrosetta_shs_%j.err
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --gres=gpu:A40:1

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
cd "$BASE_DIR"

mkdir -p logs

module purge
module load devel/cuda

CONDA_ENV="${CONDA_ENV:-trRNA}"
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

PREP_DIR="${PREP_DIR:-${BASE_DIR}/prepared_shs_inputs}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/predictions_trrosetta}"
TRROSETTA_ROOT="${TRROSETTA_ROOT:-${BASE_DIR}/../../trRosettaRNA_v1.1}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_trrosetta_shs_batch.py}"

python "$PIPELINE_SCRIPT" \
  --fasta-dir "$PREP_DIR/fasta" \
  --a3m-dir "$PREP_DIR/a3m" \
  --output-dir "$OUT_DIR" \
  --trrosetta-root "$TRROSETTA_ROOT" \
  --gpu 0 \
  --n-cpu "$((SLURM_CPUS_PER_TASK / 2))" \
  --n-decoys 5 \
  --skip-existing \
  --log-file "logs/trrosetta_shs_${SLURM_JOB_ID}.log"
