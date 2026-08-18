#!/bin/bash
#SBATCH --job-name=shs_score
#SBATCH --partition=cpu-single
#SBATCH --output=potts_shs/runs/rna_only/score_%j.out
#SBATCH --error=potts_shs/runs/rna_only/score_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb

# Scores every arm together once the per-arm fold jobs finish. Chained with
# --dependency=afterany so one failed arm does not block scoring the others; scoring is
# resumable and simply skips jobs with no DSSR annotation.
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

PY=$HOME/.conda/envs/synEvo/bin/python
WD=${WORKDIR:-potts_shs/runs/rna_only}
CFG=potts_shs/configs/config.helix.single.yaml

echo "Host: $(hostname)  JID: ${SLURM_JOB_ID:-none}  start: $(date)"

# Catch anything a per-arm job left unannotated (e.g. an arm that hit its walltime).
echo "=== annotate (sweep all arms) ==="
$PY -m potts_shs.eval.run annotate --config "$CFG" --work-dir "$WD"

echo "=== score (all pair types, sensitivity check) ==="
$PY -m potts_shs.eval.run score --config "$CFG" --work-dir "$WD"
cp "$WD/scores.csv" "$WD/scores_all_pairs.csv"

echo "=== score (canonical-only, headline) ==="
$PY -m potts_shs.eval.run score --config "$CFG" --work-dir "$WD" --canonical-only

echo "SCORE_COMPLETE  end: $(date)"
echo "  headline    : $WD/scores.csv"
echo "  sensitivity : $WD/scores_all_pairs.csv"
