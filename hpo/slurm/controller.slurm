#!/bin/bash
set -euo pipefail

af3_module="$1"
controller_python="$2"
controller_module="$3"
config_path="$4"

module load "$af3_module"

# Capture the Python selected by the AF3 module before explicitly starting the
# controller with its own environment.
export AF3_PYTHON
AF3_PYTHON="$(command -v python)"

if [[ -z "$AF3_PYTHON" ]]; then
    echo "Could not locate the AF3 Python interpreter" >&2
    exit 1
fi

exec "$controller_python" \
    -u \
    -m "$controller_module" \
    --config "$config_path"
