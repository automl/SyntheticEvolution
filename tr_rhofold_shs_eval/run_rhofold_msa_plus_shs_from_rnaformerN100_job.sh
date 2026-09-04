#!/bin/bash
#SBATCH --job-name=rhofold_msa_plus_shs_from_rnaformerN100
#SBATCH --partition=gpu-single
#SBATCH --output=logs/rhofold_msa_plus_shs_from_rnaformerN100_%j.out
#SBATCH --error=logs/rhofold_msa_plus_shs_from_rnaformerN100_%j.err
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32gb
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

PREP_DIR="${PREP_DIR:-${BASE_DIR}/prepared_msa_plus_shs_from_rnaformerN100_folder}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/predictions_rhofold_on_msa_plus_shs_from_rnaformerN100_folder}"
RHOFOLD_BIN="${RHOFOLD_BIN:-${BASE_DIR}/../../RhoFold/inference.py}"
RHOFOLD_CKPT="${RHOFOLD_CKPT:-}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_rhofold_shs_batch.py}"

CMD=(
  python "$PIPELINE_SCRIPT"
  --fasta-dir "$PREP_DIR/fasta"
  --a3m-dir "$PREP_DIR/a3m"
  --msa-mode provided
  --output-dir "$OUT_DIR"
  --rhofold-bin "$RHOFOLD_BIN"
  --device cuda
  --skip-existing
  --log-file "logs/rhofold_msa_plus_shs_from_rnaformerN100_${SLURM_JOB_ID}.log"
)

if [[ -n "$RHOFOLD_CKPT" ]]; then
  CMD+=(--checkpoint "$RHOFOLD_CKPT")
fi

"${CMD[@]}"
