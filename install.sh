#!/usr/bin/env bash
set -Eeuo pipefail

# ==========================
# Config
# ==========================
ENV_FILE="env/environment.yml"
SPOT_REPO="https://github.com/jaswindersingh2/SPOT-RNA.git"
SPOT_MODELS_PRIMARY="https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz"
SPOT_MODELS_FALLBACK="https://app.nihaocloud.com/f/fbf3315a91d542c0bdc2/?dl=1"
RNABENCH_REPO="https://github.com/automl/RnaBench.git"
EXTERNAL_DIR="external_algorithms"

log()  { printf "\033[1;32m[+]\033[0m %s\n" "$*"; }
warn() { printf "\033[1;33m[!]\033[0m %s\n" "$*"; }
die()  { printf "\033[1;31m[x]\033[0m %s\n" "$*" >&2; exit 1; }

need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command '$1' not found."; }

fetch() {
  local url="$1" out="$2"
  if command -v wget >/dev/null 2>&1; then
    wget -q --show-progress --progress=bar:force:noscroll --tries=3 -O "$out" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 --retry-delay 2 -o "$out" "$url"
  else
    die "Neither wget nor curl is available; please install one."
  fi
}

# Try very hard to discover and initialize conda/mamba
try_source_any_conda_hook() {
  # Fast path: conda provides a hook that works even if not initialized
  if command -v conda >/dev/null 2>&1; then
    # Attempt to source its hook for 'conda run' stability
    local hook
    if hook="$(conda shell.bash hook 2>/dev/null)"; then eval "$hook" && return 0; fi
  fi

  # Common install roots
  local roots=(
    "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" "$HOME/miniforge3"
    "/opt/conda" "/usr/local/conda" "/usr/local/miniconda" "/usr/local/miniconda3"
  )
  for r in "${roots[@]}"; do
    [[ -f "$r/etc/profile.d/conda.sh" ]] && { # shellcheck disable=SC1090
      source "$r/etc/profile.d/conda.sh" && return 0
    }
    [[ -f "$r/etc/profile.d/mamba.sh" ]] && { # shellcheck disable=SC1090
      source "$r/etc/profile.d/mamba.sh" && return 0
    }
  done

  # Deep scan for conda.sh/mamba.sh under $HOME and common prefixes
  shopt -s nullglob globstar
  local candidates=(
    "$HOME"/**/etc/profile.d/conda.sh
    "$HOME"/**/etc/profile.d/mamba.sh
    /opt/**/etc/profile.d/conda.sh
    /opt/**/etc/profile.d/mamba.sh
    /usr/local/**/etc/profile.d/conda.sh
    /usr/local/**/etc/profile.d/mamba.sh
  )
  for f in "${candidates[@]}"; do
    # shellcheck disable=SC1090
    source "$f" && return 0
  done
  shopt -u globstar

  return 1
}

# Install micromamba into ~/.local/bin if nothing found
bootstrap_micromamba() {
  log "No conda/mamba found — bootstrapping micromamba to ~/.local/bin ..."
  need_cmd mkdir
  need_cmd tar
  need_cmd uname

  local os arch url tmpdir bin
  os="$(uname -s | tr '[:upper:]' '[:lower:]')"
  arch="$(uname -m)"
  case "$arch" in
    x86_64|amd64) arch="64";;
    aarch64|arm64) arch="aarch64";;
    *) die "Unsupported CPU architecture: $arch";;
  esac

  case "$os" in
    linux)  url="https://micro.mamba.pm/api/micromamba/linux-$arch/latest";;
    darwin) url="https://micro.mamba.pm/api/micromamba/osx-$arch/latest";;
    *) die "Unsupported OS: $os";;
  esac

  tmpdir="$(mktemp -d)"
  bin="$HOME/.local/bin"
  mkdir -p "$bin"
  fetch "$url" "$tmpdir/micromamba.tar.bz2"
  tar -xjf "$tmpdir/micromamba.tar.bz2" -C "$tmpdir" --strip-components=1 bin/micromamba
  install -m 0755 "$tmpdir/micromamba" "$bin/micromamba"
  rm -rf "$tmpdir"
  export PATH="$bin:$PATH"
  command -v micromamba >/dev/null 2>&1 || die "micromamba bootstrap failed."
  log "micromamba installed at $bin/micromamba"
}

pick_env_mgr() {
  # Prefer mamba > micromamba > conda
  if command -v mamba >/dev/null 2>&1; then
    echo "mamba"; return
  fi
  if command -v micromamba >/dev/null 2>&1; then
    echo "micromamba"; return
  fi
  if command -v conda >/dev/null 2>&1; then
    echo "conda"; return
  fi

  # Try sourcing hooks to reveal hidden installs
  if try_source_any_conda_hook; then
    if command -v mamba >/dev/null 2>&1; then echo "mamba"; return; fi
    if command -v conda >/dev/null 2>&1; then echo "conda"; return; fi
  fi

  # Final fallback: bootstrap micromamba
  bootstrap_micromamba
  echo "micromamba"
}

env_run() {
  local mgr="$1" env="$2"; shift 2
  case "$mgr" in
    micromamba|mamba|conda) "$mgr" run -n "$env" "$@";;
    *) die "Unknown env manager: $mgr";;
  esac
}

create_or_update_env() {
  local mgr="$1" env="$2" file="$3"

  # Does env exist?
  local exists=""
  set +e
  if [[ "$mgr" == "micromamba" ]]; then
    exists="$("$mgr" env list 2>/dev/null | awk '{print $1}' | grep -Fx "$env")"
  else
    exists="$("$mgr" env list 2>/dev/null | awk '{print $1}' | grep -Fx "$env")"
  fi
  set -e

  if [[ -n "$exists" ]]; then
    log "Environment '$env' exists — updating from $file …"
    if [[ "$mgr" == "micromamba" ]]; then
      # micromamba 'install -f' acts as update; --prune removes orphans
      "$mgr" install -y -n "$env" -f "$file" --prune
    else
      "$mgr" env update -n "$env" -f "$file" --prune
    fi
  else
    log "Creating environment '$env' from $file …"
    if [[ "$mgr" == "micromamba" ]]; then
      "$mgr" create -y -n "$env" -f "$file"
    else
      "$mgr" env create -n "$env" -f "$file"
    fi
  fi
}

# ==========================
# Pre-flight
# ==========================
need_cmd git
need_cmd tar
[[ -f "$ENV_FILE" ]] || die "Environment file '$ENV_FILE' not found. Run from repo root or update ENV_FILE."

# Resolve project root (repo top if inside a git checkout)
PROJECT_ROOT="$(pwd)"
set +e
GIT_TOP="$(git rev-parse --show-toplevel 2>/dev/null)"
set -e
[[ -n "${GIT_TOP:-}" ]] && PROJECT_ROOT="$GIT_TOP"
log "Project root: $PROJECT_ROOT"

# Pick/create a manager (auto-discovery + bootstrap if needed)
ENV_MGR="$(pick_env_mgr)"
log "Using environment manager: $ENV_MGR"

# Get env name from YAML (fallback to rnabench-env)
if command -v yq >/dev/null 2>&1; then
  ENV_NAME="$(yq -r '.name // empty' "$ENV_FILE" || true)"
else
  ENV_NAME="$(awk '/^[[:space:]]*name:[[:space:]]*/ {print $2; exit}' "$ENV_FILE" || true)"
fi
[[ -n "${ENV_NAME:-}" ]] || ENV_NAME="rnabench-env"
log "Environment name: $ENV_NAME"

# Create/update env
create_or_update_env "$ENV_MGR" "$ENV_NAME" "$ENV_FILE"

# pip install -e . inside env
log "Installing current repo into '$ENV_NAME' (editable)…"
env_run "$ENV_MGR" "$ENV_NAME" python -m pip install --upgrade pip
env_run "$ENV_MGR" "$ENV_NAME" python -m pip install --no-cache-dir -e "$PROJECT_ROOT"

# External algorithms
mkdir -p "$PROJECT_ROOT/$EXTERNAL_DIR"
cd "$PROJECT_ROOT/$EXTERNAL_DIR"

# --- SPOT-RNA ---
log "Installing SPOT-RNA…"
if [[ -d SPOT-RNA/.git ]]; then
  log "SPOT-RNA already cloned — pulling latest…"
  (cd SPOT-RNA && git fetch --all --tags --prune && git pull --ff-only)
else
  git clone "$SPOT_REPO" SPOT-RNA
fi

# Models (skip if already unpacked)
if [[ ! -d SPOT-RNA/models ]]; then
  log "Fetching SPOT-RNA models…"
  cd SPOT-RNA
  TMP_TAR="SPOT-RNA-models.tar.gz"
  set +e
  if ! fetch "$SPOT_MODELS_PRIMARY" "$TMP_TAR"; then
    warn "Primary URL failed — trying fallback…"
    if ! fetch "$SPOT_MODELS_FALLBACK" "$TMP_TAR"; then
      die "Failed to download SPOT-RNA models from both URLs."
    fi
  fi
  set -e
  tar -xvzf "$TMP_TAR"
  rm -f "$TMP_TAR"
  cd ..
else
  log "SPOT-RNA models already present — skipping."
fi

# --- RnaBench ---
log "Cloning RnaBench…"
if [[ -d RnaBench/.git ]]; then
  log "RnaBench already cloned — pulling latest…"
  (cd RnaBench && git fetch --all --tags --prune && git pull --ff-only)
else
  git clone "$RNABENCH_REPO" RnaBench
fi

cd ..

mv external_algorithms/RnaBench/RnaBench .
rm -rf external_algorithms/RnaBench

log "All set!"
echo ""
echo "Use the environment via:"
echo "  $ENV_MGR run -n $ENV_NAME python -c 'import sys; print(sys.executable)'"
if [[ "$ENV_MGR" == "micromamba" ]]; then
  echo "Or activate: micromamba activate $ENV_NAME"
else
  echo "Or activate: conda activate $ENV_NAME"
fi
