#!/usr/bin/env bash
set -Eeuo pipefail

# --- resolve script directory (follow symlinks) ---
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ "$SOURCE" != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"

# --- try git root; else climb up to find sitecustomize.py; else use script dir ---
REPO_ROOT="$SCRIPT_DIR"
if command -v git >/dev/null 2>&1; then
  if GIT_TOP="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null)"; then
    REPO_ROOT="$GIT_TOP"
  fi
fi
if [[ ! -f "$REPO_ROOT/sitecustomize.py" ]]; then
  # walk up a few levels to find sitecustomize.py
  CAND="$SCRIPT_DIR"
  for _ in {1..5}; do
    [[ -f "$CAND/sitecustomize.py" ]] && { REPO_ROOT="$CAND"; break; }
    PARENT="$(dirname "$CAND")"
    [[ "$PARENT" == "$CAND" ]] && break
    CAND="$PARENT"
  done
fi

# --- choose python from current env ---
PYTHON_BIN="${PYTHON:-$(command -v python || true)}"
if [[ -z "${PYTHON_BIN}" ]]; then
  PYTHON_BIN="$(command -v python3 || true)"
fi
if [[ -z "${PYTHON_BIN}" ]]; then
  echo "[-] Could not find a Python interpreter in PATH." >&2
  exit 1
fi

# --- force CPU everywhere ---
export CUDA_VISIBLE_DEVICES=""
export HIP_VISIBLE_DEVICES=""
export PYTORCH_ENABLE_MPS_FALLBACK=1

# --- prepend repo root for sitecustomize.py ---
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

# --- move into the repo so relative paths work ---
cd "$REPO_ROOT"

# helpful heads-up if sitecustomize.py is missing
if [[ ! -f "${REPO_ROOT}/sitecustomize.py" ]]; then
  echo "[!] ${REPO_ROOT}/sitecustomize.py not found — .cuda() calls may still try GPU if visible." >&2
fi

# --- arg handling: support both 'script.py ...' and 'python script.py ...' ---
if [[ $# -eq 0 ]]; then
  echo "Usage:"
  echo "  $0 your_script.py [args...]"
  echo "  $0 python your_script.py [args...]   # also works"
  echo "  $0 -m your.module [args...]          # python -m style"
  exit 2
fi

# If first arg is literally 'python' or 'python3', drop it to avoid 'python python ...'
if [[ "$1" == "python" || "$1" == "python3" ]]; then
  shift
fi

# If first arg looks like a Python entry, run via python; otherwise exec as-is
if [[ "$1" == "-m" || "$1" == *.py ]]; then
  exec "$PYTHON_BIN" "$@"
else
  exec "$@"
fi
