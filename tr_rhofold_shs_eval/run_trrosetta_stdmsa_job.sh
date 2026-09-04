#!/bin/bash
#SBATCH --job-name=trrosetta_stdmsa
#SBATCH --partition=gpu-single
#SBATCH --output=logs/trrosetta_stdmsa_%j.out
#SBATCH --error=logs/trrosetta_stdmsa_%j.err
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

PREP_SHS_DIR="${PREP_SHS_DIR:-${BASE_DIR}/prepared_shs_inputs}"
PREP_STD_DIR="${PREP_STD_DIR:-${BASE_DIR}/prepared_standard_msa}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/predictions_trrosetta_stdmsa}"
TRROSETTA_ROOT="${TRROSETTA_ROOT:-${BASE_DIR}/../../trRosettaRNA_v1.1}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_trrosetta_shs_batch.py}"

python "$PIPELINE_SCRIPT" \
  --fasta-dir "$PREP_SHS_DIR/fasta" \
  --a3m-dir "$PREP_STD_DIR/a3m" \
  --msa-mode provided \
  --output-dir "$OUT_DIR" \
  --trrosetta-root "$TRROSETTA_ROOT" \
  --gpu 0 \
  --n-cpu "$((SLURM_CPUS_PER_TASK / 2))" \
  --n-decoys 5 \
  --skip-existing \
  --log-file "logs/trrosetta_stdmsa_${SLURM_JOB_ID}.log"
