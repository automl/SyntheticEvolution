#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <pred_data_dir> <gt_data_dir> <eval type; either rna_rna or protein_rna>"
  exit 1
}

if [[ $# -ne 3 ]]; then
  usage
fi

PRED_DIR="$1"
GT_DIR="$2"
TYPE="$3"

if [[ "$TYPE" != "rna_rna" && "$TYPE" != "protein_rna" ]]; then
  echo "Error: eval type must be 'rna_rna' or 'protein_rna', got '$TYPE'."
  usage
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$SCRIPT_DIR"
cd "$PROJECT_ROOT"

python prepare_evaluation_data.py \
  --pred_data_dir "$PRED_DIR" \
  --gt_data_dir   "$GT_DIR"

DB_FILE="database/rbpDatabase.db"
if [[ -f "$DB_FILE" ]]; then
  echo "Removing old database at $DB_FILE"
  rm "$DB_FILE"
fi

PARENT_FOLDER="data/$TYPE"
PRED_FOLDER="$PARENT_FOLDER/af"
REF_FOLDER="$PARENT_FOLDER/pdb"
PRED_TABLE="pred_$TYPE"
REF_TABLE="exp_$TYPE"
DB_PATH="$PROJECT_ROOT/$DB_FILE"

CONFIG_JSON="database/config.json"
mkdir -p "$(dirname "$CONFIG_JSON")"
cat > "$CONFIG_JSON" <<EOF
{
  "parent_folder": "$PROJECT_ROOT/$PARENT_FOLDER",
  "pred_folder":   "$PROJECT_ROOT/$PRED_FOLDER",
  "ref_folder":    "$PROJECT_ROOT/$REF_FOLDER",
  "pred_table":    "$PRED_TABLE",
  "ref_table":     "$REF_TABLE",
  "db_path":       "$DB_PATH"
}
EOF

echo "Wrote config to $CONFIG_JSON"

cp $CONFIG_JSON "$PARENT_FOLDER/config.json"

python -m evaluation.mainExecution ./data -dssr