#!/bin/bash
#SBATCH --job-name=rna_verify
#SBATCH --partition=gpu-single
#SBATCH --output=potts_shs/runs/rna_only/verify_%j.out
#SBATCH --error=potts_shs/runs/rna_only/verify_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --gres=gpu:A100:1
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
module purge; module load system/singularity
echo "start: $(date)"
singularity exec --nv --bind /gpfs/bwfor/work/ws/fr_ds595-shs-paper:/gpfs/bwfor/work/ws/fr_ds595-shs-paper --bind /gpfs/bwfor/work/ws/fr_tr1040-dominik_share:/gpfs/bwfor/work/ws/fr_tr1040-dominik_share --bind /gpfs/bwfor/work/ws/fr_ds595-shs-paper/SyntheticEvolution/potts_shs/runs/rna_only/inputs:/gpfs/bwfor/work/ws/fr_ds595-shs-paper/SyntheticEvolution/potts_shs/runs/rna_only/inputs --bind /gpfs/bwfor/work/ws/fr_ds595-shs-paper/SyntheticEvolution/potts_shs/runs/rna_only/af3:/gpfs/bwfor/work/ws/fr_ds595-shs-paper/SyntheticEvolution/potts_shs/runs/rna_only/af3 /gpfs/bwfor/work/ws/fr_tr1040-dominik_share/alphafold3.sif python /gpfs/bwfor/work/ws/fr_tr1040-dominik_share/alphafold3/run_alphafold.py --json_path=/gpfs/bwfor/work/ws/fr_ds595-shs-paper/SyntheticEvolution/potts_shs/runs/rna_only/inputs/1Q9A_nomsa_s1_a1.json --output_dir=/gpfs/bwfor/work/ws/fr_ds595-shs-paper/SyntheticEvolution/potts_shs/runs/rna_only/af3 --model_dir=/gpfs/bwfor/work/ws/fr_ds595-shs-paper/alphafold3/models --norun_data_pipeline --run_inference=True
echo "exit=$?"
N=$(find potts_shs/runs/rna_only/af3/1q9a_nomsa_s1_a1 -name '*_model.cif' 2>/dev/null | wc -l)
echo "VERIFY: model cifs produced = $N"
[ "$N" -ge 1 ] && echo "VERIFY_OK" || { echo "VERIFY_FAILED"; exit 1; }
