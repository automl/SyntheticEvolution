#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
set -euo pipefail
# submit_controller.py supplies absolute interpreter, entry point and config paths.
exec "$1" -u "$2" --config "$3"
