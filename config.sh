# Site configuration; every script sources it. Defaults need nothing but this repository and the two unzipped
# download bundles in ./downloads/. Override any variable in the environment or in config.local.sh (not tracked).

export CRC_ROOT="${CRC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$CRC_ROOT/config.local.sh" ] && source "$CRC_ROOT/config.local.sh"

# Download bundles (unzip 2026_neurips_shs_inputs.zip and 2026_neurips_shs_outputs.zip into CRC_DATA)
export CRC_DATA=${CRC_DATA:-$CRC_ROOT/downloads}
export CRC_INPUTS=${CRC_INPUTS:-$CRC_DATA/2026_neurips_shs_inputs}
export CRC_OUTPUTS=${CRC_OUTPUTS:-$CRC_DATA/2026_neurips_shs_outputs}

# Reference structures (inputs bundle, ground_truth/)
export GT_SINGLE=${GT_SINGLE:-$CRC_INPUTS/ground_truth/single_91}     # 91 mmCIF
export GT_COMPLEX=${GT_COMPLEX:-$CRC_INPUTS/ground_truth/complex_377} # 377 mmCIF
export GT_LONG=${GT_LONG:-$CRC_INPUTS/ground_truth/long_25}           # 25 mmCIF

# Metric databases: one SQLite DB per <tag>_<idset>_s<seed>_m<sample>. Default: your own re-scored DBs in results/db*
# if present, else the shipped ones from the outputs bundle (metric_dbs/benchmark, metric_dbs/long_rna).
_db() { if [ -d "$CRC_ROOT/results/$1" ] || [ ! -d "$CRC_OUTPUTS/metric_dbs/$2" ]; then echo "$CRC_ROOT/results/$1"; else echo "$CRC_OUTPUTS/metric_dbs/$2"; fi; }
export DB_DIR=${DB_DIR:-$(_db db benchmark)}
export DB_LONG=${DB_LONG:-$(_db db_long long_rna)}
export DB_RHOFOLD=${DB_RHOFOLD:-$(_db db_rhofold rhofold)}   # RhoFold+ runs (one deterministic run per input)
export WORK=${WORK:-$CRC_ROOT/work}              # eval scratch + logs (gitignored)

# Only for re-generating / re-folding: where new AF3 inputs and outputs are written. (The authors' original
# layout used dataset-named dirs below these roots; analysis code reads the bundles whenever they are present.)
export AF_INPUT=${AF_INPUT:-$WORK/af_input}
export AF_OUTPUT=${AF_OUTPUT:-$WORK/af_output}

# AlphaFold 3 v3.0.0 (NOT distributed: obtain weights from Google DeepMind, build the container yourself)
export AF3_HOME=${AF3_HOME:-$HOME/af3}           # contains alphafold3/ (v3.0.0 @ 2e56235)
export AF3_SIF=${AF3_SIF:-$AF3_HOME/alphafold3.sif}
export AF3_MODELS=${AF3_MODELS:-$AF3_HOME/models}  # af3.bin
export AF3_DB=${AF3_DB:-$AF3_HOME/public_databases} # unused at inference (inputs carry MSAs/templates), bound for completeness

# Conda environments (specs in env/)
export EVAL_ENV=${EVAL_ENV:-synEvo}        # run_eval + analysis (python 3.9, pymol-open-source 2.5)
export OST_ENV=${OST_ENV:-ost}             # OpenStructure 2.11.1 (resolved by name via conda)
export RNABENCH_ENV=${RNABENCH_ENV:-RnaBench}  # SHS generator, cd-hit-est-2d
export RNABENCH_HOME=${RNABENCH_HOME:-$HOME/RnaBench}
if [ -z "${EVAL_PY:-}" ]; then   # python of the eval env (falls back to the active python3)
  for p in "${CONDA_PREFIX:-/nonexistent}/../$EVAL_ENV/bin/python" "$HOME/.conda/envs/$EVAL_ENV/bin/python"; do
    [ -x "$p" ] && EVAL_PY=$p && break
  done
fi
export EVAL_PY=${EVAL_PY:-python3}

# Optional: Rfam CMs for RNA-family annotation (not used by any table; absent => family "None")
export RFAM_DIR=${RFAM_DIR:-$CRC_ROOT/eval/engine/evaluation/data/rfam}
