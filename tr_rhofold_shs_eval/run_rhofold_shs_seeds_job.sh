#!/bin/bash
#SBATCH --job-name=rhofold_shs_seeds
#SBATCH --partition=gpu-single
#SBATCH --output=logs/rhofold_shs_seeds_%A_%a.out
#SBATCH --error=logs/rhofold_shs_seeds_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32gb
#SBATCH --gres=gpu:A100:1
#SBATCH --array=0-2

# One array task per SHS seed (0, 42, 137): regenerate VARIANT's SHS MSAs under that
# seed, then fold them with RhoFold. RhoFold inference is deterministic, so the seed
# is applied to the SHS generator, i.e. where it actually changes anything.
#
#   sbatch --export=ALL,VARIANT=rnafold run_rhofold_shs_seeds_job.sh
#   sbatch --export=ALL,VARIANT=rnafold --array=1 run_rhofold_shs_seeds_job.sh  # seed 42 only
#   bash submit_rhofold_shs_seed_jobs.sh                                        # all variants
#
# Variants: see `python generate_shs_seed_inputs.py --list-variants`.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
cd "$BASE_DIR"

mkdir -p logs

VARIANT="${VARIANT:-}"
if [[ -z "$VARIANT" ]]; then
  echo "VARIANT is unset. Submit with e.g.:" >&2
  echo "  sbatch --export=ALL,VARIANT=rnafold $(basename "$0")" >&2
  exit 2
fi

read -r -a SEEDS <<< "${SEEDS:-0 42 137}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX >= ${#SEEDS[@]} )); then
  echo "Array index $IDX has no seed in (${SEEDS[*]})." >&2
  exit 2
fi
SEED="${SEEDS[$IDX]}"

module purge
module load devel/cuda

CONDA_ENV="${CONDA_ENV:-rhofold}"
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# The SHS generator needs RnaBench, which only loads in synEvo; call that
# interpreter directly instead of switching the whole job's env.
SYNEVO_PY="${SYNEVO_PY:-$HOME/.conda/envs/synEvo/bin/python}"
if [[ ! -x "$SYNEVO_PY" ]]; then
  echo "synEvo interpreter not found: $SYNEVO_PY (override with SYNEVO_PY=...)" >&2
  exit 2
fi
GEN_SCRIPT="${GEN_SCRIPT:-${BASE_DIR}/generate_shs_seed_inputs.py}"
SHS_GENERATOR="${SHS_GENERATOR:-${BASE_DIR}/../rna_msa_generator_base_pair.py}"
N_MSA="${N_MSA:-100}"
SHS_PARAMS="${SHS_PARAMS:-}"          # JSON dict of extra generator knobs; empty => defaults

RHOFOLD_BIN="${RHOFOLD_BIN:-${BASE_DIR}/../../RhoFold/inference.py}"
RHOFOLD_CKPT="${RHOFOLD_CKPT:-}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-${BASE_DIR}/run_rhofold_shs_batch.py}"

RUN_PREP="${RUN_PREP:-1}"
RUN_FOLD="${RUN_FOLD:-1}"

# Dir names live in generate_shs_seed_inputs.py so prep and fold cannot drift apart.
eval "$("$SYNEVO_PY" "$GEN_SCRIPT" --variant "$VARIANT" --seed "$SEED" --print-paths)"
PREP_DIR="${PREP_DIR:-$DEFAULT_PREP_DIR}"
OUT_DIR="${OUT_DIR:-$DEFAULT_OUT_DIR}"

echo "Variant : $VARIANT   (pairs method: $METHOD)"
echo "Seed    : $SEED      (array index $IDX of ${SEEDS[*]})"
echo "Inputs  : $PREP_DIR"
echo "Output  : $OUT_DIR"

if [[ "$RUN_PREP" == "1" ]]; then
  echo "--- generating SHS inputs (seed $SEED) ---"
  GEN_CMD=(
    "$SYNEVO_PY" "$GEN_SCRIPT"
    --variant "$VARIANT"
    --seed "$SEED"
    --n-msa "$N_MSA"
    --out-dir "$PREP_DIR"
    --generator "$SHS_GENERATOR"
    --gen-python "$SYNEVO_PY"
    --skip-existing
  )
  if [[ -n "$SHS_PARAMS" ]]; then
    GEN_CMD+=(--shs-params "$SHS_PARAMS")
  fi
  # A single target failing to generate should not cost the whole variant its fold,
  # so only abort when nothing at all was written.
  gen_rc=0
  "${GEN_CMD[@]}" || gen_rc=$?
  if (( gen_rc != 0 )); then
    n_a3m=$(find "$PREP_DIR/a3m" -name '*.a3m' 2>/dev/null | wc -l)
    if (( n_a3m == 0 )); then
      echo "SHS generation failed (rc=$gen_rc) and wrote no A3M; aborting." >&2
      exit "$gen_rc"
    fi
    echo "WARNING: SHS generation rc=$gen_rc; folding the $n_a3m A3Ms it did write."
  fi
fi

if [[ "$RUN_FOLD" == "1" ]]; then
  echo "--- folding with RhoFold ---"
  CMD=(
    python "$PIPELINE_SCRIPT"
    --fasta-dir "$PREP_DIR/fasta"
    --a3m-dir "$PREP_DIR/a3m"
    --msa-mode provided
    --output-dir "$OUT_DIR"
    --rhofold-bin "$RHOFOLD_BIN"
    --device cuda
    --skip-existing
    --log-file "logs/rhofold_shs_${VARIANT}_shsseed${SEED}_${SLURM_JOB_ID:-local}.log"
  )
  if [[ -n "$RHOFOLD_CKPT" ]]; then
    CMD+=(--checkpoint "$RHOFOLD_CKPT")
  fi
  "${CMD[@]}"
fi

echo "Done: variant=$VARIANT seed=$SEED"
