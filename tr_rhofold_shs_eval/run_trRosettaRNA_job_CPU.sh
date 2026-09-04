#!/bin/bash
# =============================================================================
#  trRosettaRNA prediction pipeline — SLURM job script for Helix
#
#  Prerequisites (one-time setup):
#    1. Download & unzip:
#         wget https://zenodo.org/records/10616895/files/trRosettaRNA_v1.1.zip
#         unzip trRosettaRNA_v1.1.zip
#    2. Create conda env from the bundled spec:
#         conda env create -f trRosettaRNA_v1.1/linux.yml -n trrosettarna
#         conda activate trrosettarna
#         pip install einops          # required but missing from linux.yml
#    3. Ensure PyRosetta is available in the env (needed for 3D folding step).
#       Academic licence: https://www.pyrosetta.org/downloads
#
#  Submit with:
#    sbatch trrosettarna_job.sh
#
#  Override paths at submission time:
#    sbatch --export=INPUT_DIR=/path/to/fastas,OUTPUT_DIR=/path/to/out trrosettarna_job.sh
# =============================================================================

# ── Job metadata ──────────────────────────────────────────────────────────────
#SBATCH --job-name=trrosettarna_pred
#SBATCH --partition=cpu-single         # shared CPU partition on Helix (1 node max)
#SBATCH --output=logs/trrosettarna_%j.out
#SBATCH --error=logs/trrosettarna_%j.err
#SBATCH --time=10:00:00                 # geometry step is fast; Rosetta folding is slow
#SBATCH --export=NONE                   # recommended on Helix for reproducibility

# ── Resources ─────────────────────────────────────────────────────────────────
#SBATCH --ntasks=64
#SBATCH --mem=64G

# ── Optional: email notifications ─────────────────────────────────────────────
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=you@example.org

# =============================================================================
#  USER CONFIGURATION — edit these paths before submitting
# =============================================================================

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"

# Directory containing input FASTA files (.fa or .fasta)
INPUT_DIR="${INPUT_DIR:-${BASE_DIR}/data}"

# Root output directory (one sub-folder per sequence will be created)
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/trrosettarna_predictions}"

# Full path to the Python pipeline wrapper (run_trRosettaRNA.py)
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_trRosettaRNA.py}"

# Root directory of the unzipped trRosettaRNA_v1.1 package
TRROSETTA_ROOT="${TRROSETTA_ROOT:-${BASE_DIR}/../trRosettaRNA_v1.1}"

# Conda environment name
CONDA_ENV="${CONDA_ENV:-trRNA}"

# Number of parallel Rosetta folding processes (keep ≤ cpus-per-task)
N_PROC="${N_PROC:-8}"

# =============================================================================
#  Sanity checks on package layout
# =============================================================================

PREDICT_SCRIPT="${TRROSETTA_ROOT}/predict.py"
PARAMS_DIR="${TRROSETTA_ROOT}/params"

if [[ ! -f "$PREDICT_SCRIPT" ]]; then
    echo "ERROR: trRosettaRNA predict.py not found at: $PREDICT_SCRIPT" >&2
    echo "       Make sure TRROSETTA_ROOT points to the unzipped v1.1 folder." >&2
    exit 1
fi

if [[ ! -d "$PARAMS_DIR" ]]; then
    echo "ERROR: params/ directory not found under $TRROSETTA_ROOT" >&2
    exit 1
fi

# =============================================================================
#  Environment setup
# =============================================================================

echo "======================================================"
echo "  Job ID        : $SLURM_JOB_ID"
echo "  Job name      : $SLURM_JOB_NAME"
echo "  Node          : $SLURM_NODELIST"
echo "  Submit dir    : ${SLURM_SUBMIT_DIR:-<not set>}"
echo "  CPUs          : $SLURM_CPUS_PER_TASK"
echo "  Start time    : $(date)"
echo "  Input dir     : $INPUT_DIR"
echo "  Output dir    : $OUTPUT_DIR"
echo "  trRosetta root: $TRROSETTA_ROOT"
echo "======================================================"

mkdir -p logs

# ── Load system modules ───────────────────────────────────────────────────────
module purge
module load devel/cuda                  # required on Helix before GPU programs
# module load devel/miniconda           # uncomment if conda itself needs loading

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"

# ── Activate conda environment ────────────────────────────────────────────────
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to activate conda environment '${CONDA_ENV}'" >&2
    exit 1
fi

echo "Python  : $(which python) — $(python --version)"
echo "Device  : $(python -c 'import torch; print(torch.cuda.get_device_name(0))' 2>/dev/null || echo 'GPU info unavailable')"

# =============================================================================
#  Validate inputs
# =============================================================================

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: INPUT_DIR does not exist: $INPUT_DIR" >&2
    exit 1
fi

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
    echo "ERROR: Pipeline script not found: $PIPELINE_SCRIPT" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# =============================================================================
#  Build and run the prediction command
# =============================================================================

CMD=(
    python "$PIPELINE_SCRIPT"
    --input_dir       "$INPUT_DIR"
    --output_dir      "$OUTPUT_DIR"
    --trrosetta_root  "$TRROSETTA_ROOT"
    --gpu             -1                  # GPU 0; use -1 for CPU-only; Determines which GPU is used for geometry prediction step. 
    --n_cpu           "$((${SLURM_CPUS_PER_TASK:-64} / 2))"
    --n_decoys        5
    --skip_existing                      # restart-safe
    --log_file        "logs/trrosettarna_${SLURM_JOB_ID}.log"
)

echo ""
echo "Command: ${CMD[*]}"
echo ""

"${CMD[@]}"
EXIT_CODE=$?

# =============================================================================
#  Post-run summary
# =============================================================================

echo ""
echo "======================================================"
echo "  End time   : $(date)"
echo "  Exit code  : $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
    echo "  Status     : SUCCESS"
    N_PDBS=$(find "$OUTPUT_DIR" -name "*.pdb" | wc -l)
    echo "  PDB files  : $N_PDBS"
else
    echo "  Status     : FAILED (check logs/trrosettarna_${SLURM_JOB_ID}.err)"
fi

echo "======================================================"

exit $EXIT_CODE