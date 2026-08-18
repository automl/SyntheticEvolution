#!/bin/bash
#SBATCH --job-name=shs_rna_setup
#SBATCH --partition=cpu-single
#SBATCH --output=potts_shs/runs/rna_only/setup_%j.out
#SBATCH --error=potts_shs/runs/rna_only/setup_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb

# RNA-only evaluation: the 91 entries in data/datafiles/alphafold3 that carry no protein chain.
# Plans all five arms, annotates the ground truth, then submits one GPU job per arm plus a
# scorer chained on all of them.
#
# Pairs come from data/rnaformer_predictions.pkl (per-chain; identical homodimer copies are
# deduplicated to one row per (target, sequence), giving exactly 91 targets).
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

PY=$HOME/.conda/envs/synEvo/bin/python
WD=potts_shs/runs/rna_only
CFG=potts_shs/configs/config.helix.single.yaml
IDS=$WD/rna_only_targets.txt

# Calibrated on 164 Rfam seed families at depth 100; applied to every arm so mutation rate is
# never confounded with the substitution model.
RATES="--mutation-rate-paired 0.345 --mutation-rate-unpaired 0.366 -N 100"

# Two SHS seeds: the target set is small enough to afford the replication, and it is the only
# way to tell a real difference between arms from generator sampling noise.
common=(--config "$CFG" --work-dir "$WD" --targets-pickle data/rnaformer_predictions.pkl \
        --target-ids "$IDS" --input-json-dir data/datafiles/alphafold3 --shs-seeds 1 2)

echo "Host: $(hostname)  JID: ${SLURM_JOB_ID:-none}  start: $(date)"

echo "=== arm 1/5: potts_g225 (Rfam-calibrated coupling) ==="
$PY -m potts_shs.eval.run plan "${common[@]}" --approaches potts --arm-label potts_g225 \
    --generator-flags --potts-coupling 2.25 --potts-wobble 0.40 $RATES

echo "=== arm 2/5: potts_g6 (strong coupling) ==="
$PY -m potts_shs.eval.run plan "${common[@]}" --approaches potts --arm-label potts_g6 \
    --generator-flags --potts-coupling 6.0 --potts-wobble 0.85 $RATES

echo "=== arm 3/5: watson_crick_cov (incumbent) ==="
$PY -m potts_shs.eval.run plan "${common[@]}" --approaches watson_crick_cov --generator-flags $RATES

echo "=== arm 4/5: none (depth-matched, no covariation) ==="
$PY -m potts_shs.eval.run plan "${common[@]}" --approaches none --generator-flags $RATES

# -N 1 emits the query row alone: no alignment at all, AF3 folds from its own prior.
echo "=== arm 5/5: nomsa (depth-1) ==="
$PY -m potts_shs.eval.run plan "${common[@]}" --approaches none --arm-label nomsa \
    --generator-flags --mutation-rate-paired 0.345 --mutation-rate-unpaired 0.366 -N 1

echo "=== ground-truth DSSR ==="
$PY -m potts_shs.eval.run annotate --config "$CFG" --work-dir "$WD" \
    --gt-cif-dir evaluation/predictions/gt

echo "=== submit one GPU job per arm ==="
ARM_IDS=""
for ARM in potts_g225 potts_g6 watson_crick_cov none nomsa; do
    J=$(sbatch --parsable --export=ALL,ARM=$ARM,WORKDIR=$WD --job-name=rna_$ARM \
        potts_shs/runs/rna_only/run_arm.sbatch)
    echo "  $ARM -> job $J"
    ARM_IDS="${ARM_IDS:+$ARM_IDS:}$J"
done

S=$(sbatch --parsable --dependency=afterany:"$ARM_IDS" --export=ALL,WORKDIR=$WD \
    potts_shs/runs/rna_only/score_all.sbatch)
echo "  scorer -> job $S (afterany:$ARM_IDS)"
echo "SETUP_COMPLETE  end: $(date)"
