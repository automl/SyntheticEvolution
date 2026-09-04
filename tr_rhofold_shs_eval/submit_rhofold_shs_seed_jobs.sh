#!/bin/bash
# Submits run_rhofold_shs_seeds_job.sh once per SHS variant, each a seed array
# (0/42/137 by default). Run from the login node: bash submit_rhofold_shs_seed_jobs.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-python3}"
ARRAY="${ARRAY:-0-2}"          # index into SEEDS inside the job (default 3 seeds)
DRY_RUN="${DRY_RUN:-0}"
EXTRA_EXPORT="${EXTRA_EXPORT:-}"   # e.g. "RUN_FOLD=0" to only generate the SHS inputs

VARIANTS="${VARIANTS:-$("$PYTHON" generate_shs_seed_inputs.py --list-variants | tail -n +2 | awk '{print $1}')}"

for variant in $VARIANTS; do
  export_list="ALL,VARIANT=${variant}"
  if [[ -n "$EXTRA_EXPORT" ]]; then
    export_list="${export_list},${EXTRA_EXPORT}"
  fi
  cmd=(sbatch --array="$ARRAY" --export="$export_list"
       --job-name="rhofold_shs_seeds_${variant}" run_rhofold_shs_seeds_job.sh)
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "${cmd[@]}"
  else
    "${cmd[@]}"
  fi
done
