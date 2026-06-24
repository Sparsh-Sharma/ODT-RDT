#!/usr/bin/env bash
#
#  build_odt.sh — bootstrap + build BYUignite/ODT from scratch on a clean Ubuntu.
#
#  Design principle: the C++ binary (odt.x) is the deliverable. Nothing about
#  the Python post-processing stack is allowed to block or abort the build.
#  The binary is compiled and verified FIRST; Python deps and checks come
#  afterward and only ever warn, never fail the run.
#
#  Order of operations:
#    1. Install Miniforge if conda is not already available.
#    2. Create/ensure a conda env 'odt' with the C++ BUILD dependencies:
#         Cantera 2.5.1, fmt, Boost, matched compiler/linker toolchain,
#         CMake >= 3.15, make.
#    3. Configure an out-of-source CMake build using the conda toolchain.
#    4. Build odt.x on 8 cores, copy into run/, verify libraries resolve.
#    5. ONLY THEN: add the Python post-processing stack (numpy<2, scipy,
#       matplotlib, pyyaml) and run import checks — all non-fatal.
#
#  Run from inside your ODT source tree (dir with CMakeLists.txt + src/),
#  or from anywhere — if no tree is found it clones BYUignite/ODT.
#
#  Usage:        bash build_odt.sh
#  Fresh env:    FRESH_ENV=1 bash build_odt.sh      # remove + recreate 'odt'
#
set -euo pipefail

# ----------------------------------------------------------------------------
#  config
# ----------------------------------------------------------------------------
ENV_NAME="odt"
CORES=8
MINIFORGE_DIR="${HOME}/miniforge3"
ODT_REPO="https://github.com/Sparsh-Sharma/ODT-RDT.git"
ODT_DIRNAME="$(basename "${ODT_REPO}" .git)"   # clone target dir, e.g. ODT-RDT
ODT_COMMIT=""                 # optional: pin a commit for reproducibility
FRESH_ENV="${FRESH_ENV:-0}"   # set to 1 to wipe and recreate the env
PYTHON_VERSION="3.10"         # only applied when CREATING a new env

# --- C++ BUILD dependencies (these gate the binary; must succeed) -----------
#   libcantera-devel=2.5.1 pinned: newer Cantera is not source-compatible.
#   cmake>=3.15 matches ODT's CMakeLists requirement.
BUILD_PKGS=(
  "libcantera-devel=2.5.1"
  "cmake>=3.15"
  "make"
  "cxx-compiler"
  "fmt"
  "boost-cpp"
)

# --- Python POST-PROCESSING dependencies (optional; never block the build) --
#   numpy<2 keeps the stack consistent with the Cantera-2.5.1-era ecosystem.
POST_PKGS=(
  "numpy<2"
  "scipy"
  "matplotlib"
  "pyyaml"
)

# ----------------------------------------------------------------------------
#  pretty output
# ----------------------------------------------------------------------------
if [[ -t 1 ]]; then
  BOLD=$'\e[1m'; DIM=$'\e[2m'; RED=$'\e[31m'; GRN=$'\e[32m'
  YLW=$'\e[33m'; BLU=$'\e[34m'; MAG=$'\e[35m'; CYN=$'\e[36m'; RST=$'\e[0m'
else
  BOLD=""; DIM=""; RED=""; GRN=""; YLW=""; BLU=""; MAG=""; CYN=""; RST=""
fi

banner() {
  echo
  echo "${MAG}${BOLD}  ┌─────────────────────────────────────────────┐${RST}"
  echo "${MAG}${BOLD}  │   O D T   ·   o n e - s t e p   b u i l d   │${RST}"
  echo "${MAG}${BOLD}  └─────────────────────────────────────────────┘${RST}"
  echo "${DIM}     turning turbulence into a binary, one eddy at a time${RST}"
  echo
}

step()  { echo; echo "${BLU}${BOLD}▸ $*${RST}"; }
info()  { echo "  ${CYN}·${RST} $*"; }
ok()    { echo "  ${GRN}✔${RST} $*"; }
warn()  { echo "  ${YLW}!${RST} $*"; }
die()   { echo; echo "  ${RED}${BOLD}ERROR: $*${RST}" >&2; exit 1; }

trap 'echo; echo "${RED}${BOLD}Build aborted.${RST} Check the message above."' ERR

# ----------------------------------------------------------------------------
#  spinner helpers
# ----------------------------------------------------------------------------
SPIN_FRAMES=( '⠋' '⠙' '⠹' '⠸' '⠼' '⠴' '⠦' '⠧' '⠇' '⠏' )

_spin_wait() {  # internal: $1=pid, $2=msg ; returns child exit status
  local pid="$1" msg="$2" i=0
  while kill -0 "$pid" 2>/dev/null; do
    printf "\r  ${YLW}%s${RST} %s" "${SPIN_FRAMES[i]}" "$msg"
    i=$(( (i + 1) % ${#SPIN_FRAMES[@]} ))
    sleep 0.1
  done
  wait "$pid"
}

run_spin() {  # FATAL: aborts the script if the command fails
  local msg="$1"; shift
  local log; log="$(mktemp)"
  ( "$@" ) >"$log" 2>&1 &
  if _spin_wait "$!" "$msg"; then
    printf "\r  ${GRN}✔${RST} %s\n" "$msg"; rm -f "$log"
  else
    printf "\r  ${RED}✗${RST} %s\n" "$msg"
    echo "${DIM}---- last 40 lines of output ----${RST}"
    tail -n 40 "$log" || true; rm -f "$log"
    die "step failed: $msg"
  fi
}

run_spin_soft() {  # NON-FATAL: warns and continues if the command fails
  local msg="$1"; shift
  local log; log="$(mktemp)"
  ( "$@" ) >"$log" 2>&1 &
  if _spin_wait "$!" "$msg"; then
    printf "\r  ${GRN}✔${RST} %s\n" "$msg"; rm -f "$log"; return 0
  else
    printf "\r  ${YLW}!${RST} %s ${DIM}(non-fatal — continuing)${RST}\n" "$msg"
    echo "${DIM}---- last 20 lines of output ----${RST}"
    tail -n 20 "$log" || true; rm -f "$log"; return 1
  fi
}

# ----------------------------------------------------------------------------
#  apt helper (only used if git/curl are missing on a clean box)
# ----------------------------------------------------------------------------
APT_UPDATED=0
apt_install() {
  if [[ "$APT_UPDATED" -eq 0 ]]; then sudo apt-get update; APT_UPDATED=1; fi
  sudo apt-get install -y "$@"
}

# ----------------------------------------------------------------------------
banner

# ----------------------------------------------------------------------------
#  1. locate the ODT source tree
# ----------------------------------------------------------------------------
step "Locating ODT source"
# A directory only counts as a valid tree if it actually contains the sources.
is_odt_tree() { [[ -f "$1/CMakeLists.txt" && -d "$1/src" ]]; }

if is_odt_tree "$(pwd)"; then
  ODT_DIR="$(pwd)"
  ok "Using current directory: ${BOLD}${ODT_DIR}${RST}"
else
  warn "Current directory is not an ODT tree — will use/clone ${ODT_DIRNAME}/."
  command -v git >/dev/null 2>&1 || run_spin "Installing git" apt_install git

  if [[ -d "$ODT_DIRNAME" ]] && is_odt_tree "$ODT_DIRNAME"; then
    ok "Found existing valid tree: ${BOLD}${ODT_DIRNAME}/${RST}"
  elif [[ -d "$ODT_DIRNAME" ]] && [[ -z "$(ls -A "$ODT_DIRNAME" 2>/dev/null)" ]]; then
    # directory exists but is empty (e.g. a failed prior clone) — clone into it
    warn "${ODT_DIRNAME}/ exists but is empty — cloning into it."
    run_spin "Cloning ${ODT_REPO}" git clone "$ODT_REPO" "$ODT_DIRNAME"
  elif [[ -d "$ODT_DIRNAME" ]]; then
    # directory exists, is non-empty, but is NOT a valid tree — don't guess
    die "${PWD}/${ODT_DIRNAME} exists but has no CMakeLists.txt/src. \
Remove it (rm -rf ${ODT_DIRNAME}) and re-run, or cd into your real ODT tree first."
  else
    run_spin "Cloning ${ODT_REPO}" git clone "$ODT_REPO" "$ODT_DIRNAME"
  fi

  ODT_DIR="$(pwd)/${ODT_DIRNAME}"; cd "$ODT_DIR"
  if [[ -n "$ODT_COMMIT" ]]; then
    run_spin "Checking out ${ODT_COMMIT}" git checkout "$ODT_COMMIT"
  fi
  ok "Using ODT tree: ${BOLD}${ODT_DIR}${RST}"
fi
cd "$ODT_DIR"
is_odt_tree "$ODT_DIR" || die "Not a valid ODT tree (missing CMakeLists.txt or src/): ${ODT_DIR}"

# ----------------------------------------------------------------------------
#  2. ensure conda / Miniforge
# ----------------------------------------------------------------------------
step "Checking for conda"
CONDA_SH=""
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"
  ok "Found conda at ${BOLD}${CONDA_BASE}${RST}"
elif [[ -f "${MINIFORGE_DIR}/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="${MINIFORGE_DIR}/etc/profile.d/conda.sh"
  ok "Found existing Miniforge at ${BOLD}${MINIFORGE_DIR}${RST}"
else
  info "conda not found — installing Miniforge."
  ARCH="$(uname -m)"
  case "$ARCH" in
    x86_64)  MF="Miniforge3-Linux-x86_64.sh" ;;
    aarch64) MF="Miniforge3-Linux-aarch64.sh" ;;
    *)       die "Unsupported architecture: $ARCH" ;;
  esac
  command -v curl >/dev/null 2>&1 || run_spin "Installing curl" apt_install curl
  run_spin "Downloading Miniforge (${ARCH})" \
    curl -fsSL -o "/tmp/${MF}" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/${MF}"
  run_spin "Installing Miniforge to ${MINIFORGE_DIR}" \
    bash "/tmp/${MF}" -b -p "${MINIFORGE_DIR}"
  rm -f "/tmp/${MF}"
  CONDA_SH="${MINIFORGE_DIR}/etc/profile.d/conda.sh"
  ok "Miniforge installed."
fi
[[ -f "$CONDA_SH" ]] || die "Could not find conda activation script: ${CONDA_SH}"
# shellcheck disable=SC1090
source "$CONDA_SH"

# ----------------------------------------------------------------------------
#  3. create / ensure the env with C++ BUILD deps only
# ----------------------------------------------------------------------------
step "Preparing the '${ENV_NAME}' conda environment (build deps)"

ENV_EXISTS=0
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then ENV_EXISTS=1; fi

if [[ "$FRESH_ENV" == "1" && "$ENV_EXISTS" == "1" ]]; then
  info "FRESH_ENV=1 — removing existing '${ENV_NAME}' for a clean rebuild."
  run_spin "Removing old env" conda env remove -n "${ENV_NAME}" -y
  ENV_EXISTS=0
fi

if [[ "$ENV_EXISTS" == "1" ]]; then
  # Reuse as-is. Do NOT pin a python version here: pinning against an env that
  # already has a different python is either ignored (silent) or forces a
  # churny rebuild. Just make sure the build deps are present.
  ok "Env '${ENV_NAME}' exists — ensuring build dependencies are present."
  warn "Reusing existing env; pass FRESH_ENV=1 to recreate from scratch."
  run_spin "Solving + installing build deps" \
    conda install -n "${ENV_NAME}" --override-channels -c conda-forge -y "${BUILD_PKGS[@]}"
else
  info "Creating env (python=${PYTHON_VERSION}) with build deps:"
  for pkg in "${BUILD_PKGS[@]}"; do echo "    ${DIM}-${RST} ${pkg}"; done
  run_spin "Solving + downloading build deps" \
    conda create -n "${ENV_NAME}" --override-channels -c conda-forge -y \
      "python=${PYTHON_VERSION}" "${BUILD_PKGS[@]}"
fi

conda activate "${ENV_NAME}"
[[ -n "${CONDA_PREFIX:-}" ]]        || die "CONDA_PREFIX not set after activation."
[[ -n "${CC:-}" && -n "${CXX:-}" ]] || die "Conda compilers CC/CXX not set after activation."
ok "Activated env: ${BOLD}${CONDA_PREFIX}${RST}"
ok "Toolchain: ${BOLD}$(basename "$CXX")${RST}  ($("$CXX" -dumpversion))"

# ----------------------------------------------------------------------------
#  4. configure (clean, out-of-source, conda toolchain)
# ----------------------------------------------------------------------------
step "Configuring the build"
rm -rf build
run_spin "Running CMake configure" \
  cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
ok "Build files generated in ${BOLD}build/${RST}"

# ----------------------------------------------------------------------------
#  5. build  (THE deliverable — streamed live)
# ----------------------------------------------------------------------------
step "Compiling odt.x on ${BOLD}${CORES}${RST} cores"
echo "${DIM}  live compiler output below${RST}"; echo
cmake --build build -j"${CORES}"
echo
ok "Compilation finished."

# ----------------------------------------------------------------------------
#  6. locate binary, copy to run/, verify libraries
# ----------------------------------------------------------------------------
step "Finalizing executable"
BIN=""
for cand in build/src/odt.x build/odt.x; do
  [[ -x "$cand" ]] && BIN="$cand" && break
done
[[ -n "$BIN" ]] || die "Could not find built odt.x. Check the build output above."
ok "Built executable: ${BOLD}${BIN}${RST}"

if [[ -d run ]]; then
  cp -f "$BIN" run/ && ok "Copied to ${BOLD}run/odt.x${RST}"
else
  warn "No run/ directory found; binary not copied."
fi

step "Checking shared-library resolution"
if command -v ldd >/dev/null 2>&1; then
  if ldd "$BIN" | grep -qi "not found"; then
    echo "${DIM}---- ldd ----${RST}"; ldd "$BIN" || true
    die "Some libraries did not resolve. Run inside: conda activate ${ENV_NAME}"
  fi
  ldd "$BIN" | grep -Ei 'cantera|yaml|boost|fmt' | sed 's/^/  /' || true
  ok "Library check passed."
else
  warn "ldd not found; skipped library check."
fi

# ============================================================================
#  >>> Binary is built and verified. Everything below is OPTIONAL and
#      NON-FATAL: missing Python pieces only warn, they never abort. <<<
# ============================================================================

# ----------------------------------------------------------------------------
#  7. Python post-processing stack (optional, non-fatal)
# ----------------------------------------------------------------------------
step "Adding Python post-processing stack ${DIM}(optional)${RST}"
info "These are for postprocessing only and never block the build."
POST_OK=1
if ! run_spin_soft "Solving + installing post-processing deps" \
       conda install -n "${ENV_NAME}" --override-channels -c conda-forge -y "${POST_PKGS[@]}"; then
  POST_OK=0
  warn "Post-processing deps did not fully install — the binary is still fine."
  warn "Retry later: conda install -n ${ENV_NAME} -c conda-forge \"numpy<2\" scipy matplotlib pyyaml"
fi

# ----------------------------------------------------------------------------
#  8. version report + import check (optional, non-fatal)
# ----------------------------------------------------------------------------
step "Environment summary"
cmake --version | head -n1 | sed 's/^/  /' || true
make  --version | head -n1 | sed 's/^/  /' || true

if [[ "${POST_OK}" == "1" ]]; then
  if python - <<'PY' 2>/dev/null
import sys, numpy, scipy, matplotlib, yaml
print("  Python:", sys.version.split()[0])
print("  numpy:", numpy.__version__, "| scipy:", scipy.__version__,
      "| matplotlib:", matplotlib.__version__, "| pyyaml:", yaml.__version__)
assert int(numpy.__version__.split('.')[0]) < 2, "numpy must be <2 for this stack"
PY
  then
    ok "Python post-processing stack imports cleanly."
  else
    warn "Python import check failed — postprocessing may need attention (build unaffected)."
  fi
else
  warn "Skipping Python import check (post deps not installed)."
fi

# ----------------------------------------------------------------------------
#  9. environment snapshots (optional, non-fatal)
# ----------------------------------------------------------------------------
step "Writing environment snapshots"
run_spin_soft "Exporting odt-env-history.yml" \
  bash -c "conda env export -n '${ENV_NAME}' --from-history > odt-env-history.yml" || true
run_spin_soft "Exporting odt-env-full.yml" \
  bash -c "conda env export -n '${ENV_NAME}' > odt-env-full.yml" || true

# ----------------------------------------------------------------------------
#  done
# ----------------------------------------------------------------------------
echo
echo "${GRN}${BOLD}  ✔ ODT built successfully.${RST}"
echo
echo "  ${BOLD}Run it:${RST}"
echo "    ${CYN}conda activate ${ENV_NAME}${RST}     ${DIM}# binary needs the env's libs on its rpath${RST}"
echo "    ${CYN}cd ${ODT_DIR}/run${RST}"
echo "    ${CYN}./odt.x${RST}"
echo
echo "  ${DIM}Recreate this env from scratch any time with:  FRESH_ENV=1 bash build_odt.sh${RST}"
echo
