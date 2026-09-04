#!/bin/bash
# =============================================================================
#  RhoFold prediction pipeline — SLURM job script for Helix
#
#  Submit with:
#    sbatch rhofold_job.sh
#
#  Or override paths at submission time:
#    sbatch --export=INPUT_DIR=/path/to/fastas,OUTPUT_DIR=/path/to/out rhofold_job.sh
# =============================================================================
 
# ── Job metadata ──────────────────────────────────────────────────────────────
#SBATCH --job-name=rhofold_pred
#SBATCH --partition=gpu-single          # correct Helix GPU partition (shared, 1 node max)
#SBATCH --output=logs/rhofold_%j.out
#SBATCH --error=logs/rhofold_%j.err
#SBATCH --time=06:00:00                 # max 120:00:00 on gpu-single
#SBATCH --export=NONE                   # recommended on Helix for reproducibility
 
# ── Resources ─────────────────────────────────────────────────────────────────
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb                     # Helix uses lowercase 'gb'
#SBATCH --gres=gpu:A100:1              # gpu type required on gpu-single; A100 for tensor cores
                                        # alternatives: gpu:A40:1 (FP32), gpu:H200:1

# ── Optional: email notifications ─────────────────────────────────────────────
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=dom.scheuer@gmail.com

# =============================================================================
#  USER CONFIGURATION — edit these paths before submitting
# =============================================================================

# Base directory for defaults:
# - Prefer SLURM_SUBMIT_DIR (directory from which sbatch was called)
# - Fallback to script directory if not running under Slurm
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"

# Directory containing your input FASTA files
INPUT_DIR="${INPUT_DIR:-${BASE_DIR}/data}"
echo "Input directory: $INPUT_DIR"
# Root directory where predictions will be saved (one sub-folder per sequence)
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/rhofold_predictions}"

# Full path to this Python pipeline script
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_rhofold.py}"

# RhoFold command/script:
# - If RHOFOLD_BIN is exported, use it.
# - Else prefer local repo script at ${BASE_DIR}/RhoFold/inference.py.
# - Else fallback to console command 'rhofold' from active environment.
if [[ -z "${RHOFOLD_BIN:-}" ]]; then
    if [[ -f "${BASE_DIR}/../RhoFold/inference.py" ]]; then
        RHOFOLD_BIN="${BASE_DIR}/../RhoFold/inference.py"
    else
        RHOFOLD_BIN="rhofold"
    fi
fi

# If using local RhoFold repo inference script, ensure bundled BLASTN is executable.
if [[ "$RHOFOLD_BIN" == */inference.py && -f "$RHOFOLD_BIN" ]]; then
    RHOFOLD_ROOT="$(cd -- "$(dirname -- "$RHOFOLD_BIN")" && pwd)"
    LOCAL_BLASTN="${RHOFOLD_ROOT}/rhofold/data/bin/blastn"
    if [[ -f "$LOCAL_BLASTN" && ! -x "$LOCAL_BLASTN" ]]; then
        echo "Fixing execute permission on bundled blastn: $LOCAL_BLASTN"
        chmod u+x "$LOCAL_BLASTN" 2>/dev/null || chmod a+rx "$LOCAL_BLASTN" 2>/dev/null || true
    fi
    if [[ -f "$LOCAL_BLASTN" && ! -x "$LOCAL_BLASTN" ]]; then
        echo "ERROR: Bundled blastn is not executable: $LOCAL_BLASTN" >&2
        echo "Please run: chmod +x '$LOCAL_BLASTN'" >&2
        exit 1
    fi
    FORCE_SINGLE_SEQ_PRED=0
fi

# [Optional] Path to a specific RhoFold checkpoint file
# Leave empty to use RhoFold's built-in default checkpoint
RHOFOLD_CKPT="${RHOFOLD_CKPT:-}"

# Conda environment name that has RhoFold installed
CONDA_ENV="${CONDA_ENV:-rhofold}"

# =============================================================================
#  Environment setup
# =============================================================================
 
echo "======================================================"
echo "  Job ID       : $SLURM_JOB_ID"
echo "  Job name     : $SLURM_JOB_NAME"
echo "  Node         : $SLURM_NODELIST"
echo "  Submit dir   : ${SLURM_SUBMIT_DIR:-<not set>}"
echo "  CPUs         : $SLURM_CPUS_PER_TASK"
echo "  Start time   : $(date)"
echo "  Input dir    : $INPUT_DIR"
echo "  Output dir   : $OUTPUT_DIR"
echo "  RhoFold cmd  : $RHOFOLD_BIN"
echo "======================================================"
 
# Create log directory if it doesn't exist
mkdir -p logs
 
# ── Load system modules ───────────────────────────────────────────────────────
module purge
module load devel/cuda                  # required on Helix before GPU programs
# module load devel/miniconda           # uncomment if conda itself needs loading
 
# Set OpenMP threads to match requested CPUs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
 
# ── Activate conda environment ────────────────────────────────────────────────
# Adjust the path to your conda installation:
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
    --input_dir  "$INPUT_DIR"
    --output_dir "$OUTPUT_DIR"
    --rhofold_bin "$RHOFOLD_BIN"
    --device     cuda
    --skip_existing          # restart-safe: won't redo finished predictions
    --split_multi            # process each record in multi-sequence FASTAs separately
    --log_file   "logs/rhofold_${SLURM_JOB_ID}.log"
)

if [[ "${FORCE_SINGLE_SEQ_PRED:-0}" -eq 1 ]]; then
    CMD+=(--rhofold_args --single_seq_pred True)
fi
 
# Attach checkpoint if provided
if [[ -n "$RHOFOLD_CKPT" ]]; then
    CMD+=(--checkpoint "$RHOFOLD_CKPT")
fi
 
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
echo "  End time  : $(date)"
echo "  Exit code : $EXIT_CODE"
 
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "  Status    : SUCCESS"
    # Count output PDB files as a quick sanity check
    N_PDBS=$(find "$OUTPUT_DIR" -name "*.pdb" | wc -l)
    echo "  PDB files : $N_PDBS"
else
    echo "  Status    : FAILED (check logs/rhofold_${SLURM_JOB_ID}.err)"
fi
 
echo "======================================================"
 
exit $EXIT_CODE