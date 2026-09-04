#!/bin/bash
#SBATCH --job-name=rhofold_stdmsa
#SBATCH --partition=gpu-single
#SBATCH --output=logs/rhofold_stdmsa_%j.out
#SBATCH --error=logs/rhofold_stdmsa_%j.err
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --gres=gpu:A100:1

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
cd "$BASE_DIR"

mkdir -p logs

module purge
module load devel/cuda

CONDA_ENV="${CONDA_ENV:-rhofold}"
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

PREP_DIR="${PREP_DIR:-${BASE_DIR}/prepared_db_msa_inputs}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/predictions_rhofold_stdmsa_${SLURM_JOB_ID}}"
RHOFOLD_ROOT="${RHOFOLD_ROOT:-${BASE_DIR}/../../RhoFold}"
RHOFOLD_BIN="${RHOFOLD_BIN:-${RHOFOLD_ROOT}/inference.py}"
RHOFOLD_CKPT="${RHOFOLD_CKPT:-}"
DATABASE_DPATH="${DATABASE_DPATH:-${RHOFOLD_ROOT}/database}"
BINARY_DPATH="${BINARY_DPATH:-${RHOFOLD_ROOT}/rhofold/data/bin}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_rhofold_shs_batch.py}"

export RHOFOLD_BLAST_THREADS=$SLURM_CPUS_PER_TASK
echo "Running Blast with $RHOFOLD_BLAST_THREADS threads"


CMD=(
  python "$PIPELINE_SCRIPT"
  --fasta-dir "$PREP_DIR/fasta"
  --msa-mode auto_db
  --database-dpath "$DATABASE_DPATH"
  --binary-dpath "$BINARY_DPATH"
  --output-dir "$OUT_DIR"
  --rhofold-bin "$RHOFOLD_BIN"
  --device cuda
  --skip-existing
  --log-file "logs/rhofold_stdmsa_${SLURM_JOB_ID}.log"
)

if [[ -n "$RHOFOLD_CKPT" ]]; then
  CMD+=(--checkpoint "$RHOFOLD_CKPT")
fi

"${CMD[@]}"
