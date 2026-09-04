#!/bin/bash
#SBATCH --job-name=rhofold_9CXF_retry
#SBATCH --partition=gpu-single
#SBATCH --output=logs/rhofold_9CXF_retry_%j.out
#SBATCH --error=logs/rhofold_9CXF_retry_%j.err
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --gres=gpu:A40:1

# Previous run failed with OpenCL CL_OUT_OF_RESOURCES during Amber relaxation.
# Fix: force OpenMM to use CUDA platform instead of OpenCL.

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

# Force OpenMM to use CUDA backend instead of OpenCL to avoid CL_OUT_OF_RESOURCES.
export OPENMM_DEFAULT_PLATFORM=CUDA

# OpenMM CUDA JIT kernel cache defaults to an unwritable scratch path on this cluster.
# Redirect to job-local temp dir.
export CUDA_CACHE_PATH="${SLURM_TMPDIR:-${HOME}/.cache/cuda_kernels}/openmm_${SLURM_JOB_ID}"
mkdir -p "$CUDA_CACHE_PATH"

PREP_DIR="${PREP_DIR:-${BASE_DIR}/prepared_shs_inputs_from_spotrna_folder}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/predictions_rhofold_on_shs_from_spotrna_folder}"
RHOFOLD_BIN="${RHOFOLD_BIN:-${BASE_DIR}/../../RhoFold/inference.py}"
RHOFOLD_CKPT="${RHOFOLD_CKPT:-}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_rhofold_shs_batch.py}"

IDS_FILE=$(mktemp)
echo "9CXF" > "$IDS_FILE"
trap 'rm -f "$IDS_FILE"' EXIT

# No --skip-existing: unrelaxed_model.pdb already exists from the failed run
# and would otherwise cause the batch runner to skip this target.
CMD=(
  python "$PIPELINE_SCRIPT"
  --fasta-dir "$PREP_DIR/fasta"
  --a3m-dir "$PREP_DIR/a3m"
  --msa-mode provided
  --output-dir "$OUT_DIR"
  --rhofold-bin "$RHOFOLD_BIN"
  --device cuda
  --ids-file "$IDS_FILE"
  --log-file "logs/rhofold_9CXF_retry_${SLURM_JOB_ID}.log"
)

if [[ -n "$RHOFOLD_CKPT" ]]; then
  CMD+=(--checkpoint "$RHOFOLD_CKPT")
fi

"${CMD[@]}"
