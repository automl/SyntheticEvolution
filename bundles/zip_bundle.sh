#!/bin/bash
# Build a download bundle (symlink tree) and zip it (zip stores the link targets' contents).
# Usage: sbatch --export=CRC_ROOT=<crc> bundles/zip_bundle.sh <inputs|outputs> <dest_parent_dir>
#   -> <dest_parent_dir>/2026_neurips_shs_<which>/  and  <dest_parent_dir>/2026_neurips_shs_<which>.zip (+ .sha256)
#SBATCH --job-name=crc_zip
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=72:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
source "$CRC_ROOT/config.sh"
WHICH=$1; PARENT=$2; NAME=2026_neurips_shs_$WHICH
PY=$HOME/.conda/envs/$EVAL_ENV/bin/python
rm -rf "${PARENT:?}/$NAME"
PYTHONDONTWRITEBYTECODE=1 "$PY" "$CRC_ROOT/bundles/build_bundles.py" --which "$WHICH" --dest "$PARENT/$NAME"
cd "$PARENT"
rm -f "$NAME.zip"
if [ "$WHICH" = outputs ]; then
  # slim bundle: structures + summary confidences + metric DBs; the full per-sample confidences.json
  # (PAE/atom pLDDT matrices, ~80% of the bytes) and the *_data.json input copies are left out.
  zip -r -q -6 "$NAME.zip" "$NAME" -i '*/model.cif' '*/summary_confidences.json' '*/ranking_scores.csv' \
    '*/TERMS_OF_USE.md' "$NAME/rhofold/*" "$NAME/metric_dbs/*" "$NAME/MANIFEST.tsv" "$NAME/README.txt"
else
  zip -r -q -6 "$NAME.zip" "$NAME"        # follows symlinks; zip64 for >4 GiB
fi
zip -T "$NAME.zip"
sha256sum "$NAME.zip" > "$NAME.zip.sha256"
ls -l "$NAME.zip"; cat "$NAME.zip.sha256"
echo "ZIP_DONE $NAME"
