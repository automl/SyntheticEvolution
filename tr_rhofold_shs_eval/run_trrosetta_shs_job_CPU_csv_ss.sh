#!/bin/bash
#SBATCH --job-name=trrosetta_shs_csv_ss
#SBATCH --partition=cpu-single
#SBATCH --output=logs/trrosetta_shs_csv_ss_%j.out
#SBATCH --error=logs/trrosetta_shs_csv_ss_%j.err
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64gb

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
cd "$BASE_DIR"

mkdir -p logs

module purge
module load devel/cuda

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"

CONDA_ENV="${CONDA_ENV:-trRNA}"
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

PREP_DIR="${PREP_DIR:-${BASE_DIR}/prepared_shs_inputs_from_rnaformerN100_folder}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/predictions_trrosetta_on_csv_custom_ss}"
TRROSETTA_ROOT="${TRROSETTA_ROOT:-${BASE_DIR}/../../trRosettaRNA_v1.1}"
PAIRS_CSV="${PAIRS_CSV:-${BASE_DIR}/RNA_Monomers_base_pairs.csv}"
CSV_METHOD="${CSV_METHOD:-rnaformer_pred}"
PAIRS_COLUMN="${PAIRS_COLUMN:-pairs}"
SS_SUFFIX="${SS_SUFFIX:-.ssprob.txt}"
SS_FMT="${SS_FMT:-spot_prob}"
REQUIRE_SS="${REQUIRE_SS:-1}"
OVERWRITE_SS="${OVERWRITE_SS:-0}"
IDS_FILE="${IDS_FILE:-}"

CONVERT_SCRIPT="${CONVERT_SCRIPT:-${BASE_DIR}/prepare_trrosetta_ss_from_base_pairs_csv.py}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_trrosetta_shs_batch.py}"

SS_DIR_DEFAULT="${PREP_DIR}/custom_ss_${CSV_METHOD}_${PAIRS_COLUMN}"
SS_DIR="${SS_DIR:-${SS_DIR_DEFAULT}}"

# trRosettaRNA predict.py mutates files under TRROSETTA_ROOT when invoking
# SPOT-RNA. Running multiple SLURM jobs against one shared tree can race.
TRROSETTA_USE_JOB_LOCAL_COPY="${TRROSETTA_USE_JOB_LOCAL_COPY:-1}"
if [[ "${TRROSETTA_USE_JOB_LOCAL_COPY}" == "1" ]]; then
  TRROSETTA_LOCAL_ROOT_BASE="${TRROSETTA_LOCAL_ROOT_BASE:-${SLURM_TMPDIR:-/tmp/${USER}}}"
  JOB_LOCAL_TRROSETTA_ROOT="${TRROSETTA_LOCAL_ROOT_BASE}/trRosettaRNA_v1.1_job_${SLURM_JOB_ID:-$$}"
  mkdir -p "${TRROSETTA_LOCAL_ROOT_BASE}"
  if [[ ! -d "${JOB_LOCAL_TRROSETTA_ROOT}" ]]; then
    cp -a "${TRROSETTA_ROOT}" "${JOB_LOCAL_TRROSETTA_ROOT}"
  fi
  TRROSETTA_ROOT="${JOB_LOCAL_TRROSETTA_ROOT}"
fi

CONVERT_CMD=(
  python "$CONVERT_SCRIPT"
  --pairs-csv "$PAIRS_CSV"
  --fasta-dir "$PREP_DIR/fasta"
  --method "$CSV_METHOD"
  --pairs-column "$PAIRS_COLUMN"
  --output-dir "$SS_DIR"
  --output-suffix "$SS_SUFFIX"
)

if [[ -n "$IDS_FILE" ]]; then
  CONVERT_CMD+=(--ids-file "$IDS_FILE")
fi
if [[ "$OVERWRITE_SS" == "1" ]]; then
  CONVERT_CMD+=(--overwrite)
fi

"${CONVERT_CMD[@]}"

RUN_CMD=(
  python "$PIPELINE_SCRIPT"
  --fasta-dir "$PREP_DIR/fasta"
  --a3m-dir "$PREP_DIR/a3m"
  --output-dir "$OUT_DIR"
  --trrosetta-root "$TRROSETTA_ROOT"
  --gpu -1
  --n-cpu "$((SLURM_CPUS_PER_TASK / 2))"
  --n-decoys 5
  --skip-existing
  --ss-dir "$SS_DIR"
  --ss-fmt "$SS_FMT"
  --ss-suffix "$SS_SUFFIX"
  --log-file "logs/trrosetta_shs_csv_ss_${SLURM_JOB_ID}.log"
)

if [[ "$REQUIRE_SS" == "1" ]]; then
  RUN_CMD+=(--require-ss)
fi
if [[ -n "$IDS_FILE" ]]; then
  RUN_CMD+=(--ids-file "$IDS_FILE")
fi

"${RUN_CMD[@]}"
