#!/bin/bash
set -euo pipefail

# Arguments:
#   1: AF3 module
#   2: controller Python
#   3: controller script
#   4: configuration path

af3_module="$1"
controller_python="$2"
controller_script="$3"
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
