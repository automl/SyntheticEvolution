#!/bin/bash
#SBATCH --job-name=shs_arm
#SBATCH --partition=gpu-single
#SBATCH --output=potts_shs/runs/rna_only/arm_%j.out
#SBATCH --error=potts_shs/runs/rna_only/arm_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --gres=gpu:A100:1

# One arm of the evaluation, folded sequentially on a single A100. Submit with
#   sbatch --export=ALL,ARM=<arm> --job-name=shs_<arm> potts_shs/runs/rna_only/run_arm.sbatch
#
# Three of these run in parallel against one shared work-dir: --arms narrows the manifest, and
# each fold writes into its own af3/<job_id>/ directory, so the arms never touch the same files.
# 432 targets x 2 AF3 seeds = 864 folds per arm, ~24 h at the ~1.7 min/fold measured on the smoke
# run. Resumable: resubmit verbatim and it continues from the last completed fold.
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

: "${ARM:?set ARM, e.g. sbatch --export=ALL,ARM=potts_g225 ...}"

module purge
module load system/singularity

PY=$HOME/.conda/envs/synEvo/bin/python
WD=${WORKDIR:-potts_shs/runs/rna_only}
CFG=potts_shs/configs/config.helix.single.yaml   # scheduler.kind: none -> serial, in-process

echo "Host: $(hostname)  JID: ${SLURM_JOB_ID:-none}  ARM: $ARM  start: $(date)"

echo "=== fold (sequential) ==="
$PY -m potts_shs.eval.run fold --config "$CFG" --work-dir "$WD" --arms "$ARM"

# Annotate only this arm; DSSR writes one file per job, so the arms do not collide.
echo "=== annotate ==="
$PY -m potts_shs.eval.run annotate --config "$CFG" --work-dir "$WD" --arms "$ARM"

echo "ARM_COMPLETE $ARM  end: $(date)"
