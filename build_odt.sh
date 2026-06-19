#!/usr/bin/env bash
#
#  build_odt.sh — bootstrap + build BYUignite/ODT from scratch on a clean Ubuntu.
#
#  What it does, in order:
#    1. Installs Miniforge (conda) if it isn't already on the system.
#    2. Creates a self-contained conda env 'odt' with the full ODT dependency
#       stack pinned to a known-good solve (Cantera 2.5.1, fmt, Boost, a matched
#       compiler + linker toolchain, cmake, make).
#    3. Configures an out-of-source CMake build using the conda toolchain
#       (so compiler, linker, sysroot and libs are all internally consistent).
#    4. Builds odt.x on 8 cores and drops the binary into run/.
#    5. Sanity-checks that the executable's libraries resolve.
#
#  Run it from inside your ODT source tree (the directory with CMakeLists.txt),
#  or from anywhere — if no CMakeLists.txt is found it will clone the repo.
#
#  Usage:   bash build_odt.sh
#  Re-run:  safe to run again; it reuses what already exists.
#
set -euo pipefail

# ----------------------------------------------------------------------------
#  config
# ----------------------------------------------------------------------------
ENV_NAME="odt"
CORES=8
MINIFORGE_DIR="${HOME}/miniforge3"
ODT_REPO="https://github.com/BYUignite/ODT.git"
# Pinned dependency set — one solve so fmt/Boost land at versions compatible
# with Cantera 2.5.1, and the compiler/linker come as a matched pair.
CONDA_PKGS=(
  "libcantera-devel=2.5.1"
  "cmake"
  "make"
  "cxx-compiler"
  "fmt"
  "boost-cpp"
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
die()   { echo; echo "  ${RED}${BOLD}�’✗ $*${RST}" >&2; exit 1; }

# Spinner: run a quiet long command, twirl while it works, keep a log.
SPIN_FRAMES=( '⠋' '⠙' '⠹' '⠸' '⠼' '⠴' '⠦' '⠧' '⠇' '⠏' )
run_spin() {
  local msg="$1"; shift
  local log; log="$(mktemp)"
  ( "$@" ) >"$log" 2>&1 &
  local pid=$!
  local i=0
  while kill -0 "$pid" 2>/dev/null; do
    printf "\r  ${YLW}%s${RST} %s" "${SPIN_FRAMES[i]}" "$msg"
    i=$(( (i + 1) % ${#SPIN_FRAMES[@]} ))
    sleep 0.1
  done
  if wait "$pid"; then
    printf "\r  ${GRN}✔${RST} %s\n" "$msg"
    rm -f "$log"
  else
    printf "\r  ${RED}✗${RST} %s\n" "$msg"
    echo "${DIM}---- last 30 lines of output ----${RST}"
    tail -n 30 "$log" || true
    rm -f "$log"
    die "step failed: $msg"
  fi
}

trap 'echo; echo "${RED}${BOLD}Build aborted.${RST} Check the message above."' ERR

# ----------------------------------------------------------------------------
banner

# ----------------------------------------------------------------------------
#  0. locate the ODT source tree
# ----------------------------------------------------------------------------
step "Locating ODT source"
if [[ -f "CMakeLists.txt" && -d "src" ]]; then
  ODT_DIR="$(pwd)"
  ok "Using current directory: ${BOLD}${ODT_DIR}${RST}"
else
  warn "No CMakeLists.txt here — will clone a fresh copy."
  command -v git >/dev/null 2>&1 || run_spin "Installing git" sudo apt-get install -y git
  if [[ ! -d "ODT" ]]; then
    run_spin "Cloning BYUignite/ODT" git clone "$ODT_REPO" ODT
  fi
  ODT_DIR="$(pwd)/ODT"
  ok "Using cloned tree: ${BOLD}${ODT_DIR}${RST}"
fi
cd "$ODT_DIR"

# ----------------------------------------------------------------------------
#  1. ensure conda (Miniforge) exists
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
  info "conda not found — installing Miniforge (this is the only system-level step)."
  ARCH="$(uname -m)"
  case "$ARCH" in
    x86_64)  MF="Miniforge3-Linux-x86_64.sh" ;;
    aarch64) MF="Miniforge3-Linux-aarch64.sh" ;;
    *)       die "Unsupported architecture: $ARCH" ;;
  esac
  command -v curl >/dev/null 2>&1 || run_spin "Installing curl" sudo apt-get install -y curl
  run_spin "Downloading Miniforge ($ARCH)" \
    curl -fsSL -o "/tmp/${MF}" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/${MF}"
  run_spin "Installing Miniforge to ${MINIFORGE_DIR}" \
    bash "/tmp/${MF}" -b -p "${MINIFORGE_DIR}"
  rm -f "/tmp/${MF}"
  CONDA_SH="${MINIFORGE_DIR}/etc/profile.d/conda.sh"
  ok "Miniforge installed."
fi

# Make conda usable inside this non-interactive script.
# shellcheck disable=SC1090
source "$CONDA_SH"

# ----------------------------------------------------------------------------
#  2. create / verify the 'odt' environment
# ----------------------------------------------------------------------------
step "Preparing the '${ENV_NAME}' conda environment"
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  ok "Environment '${ENV_NAME}' already exists — ensuring all packages are present."
  run_spin "Solving + installing dependencies" \
    conda install -n "${ENV_NAME}" -c conda-forge -y "${CONDA_PKGS[@]}"
else
  info "Creating env with: ${CONDA_PKGS[*]}"
  run_spin "Solving + downloading dependencies (grab a coffee)" \
    conda create -n "${ENV_NAME}" -c conda-forge -y "${CONDA_PKGS[@]}"
fi

conda activate "${ENV_NAME}"
[[ -n "${CC:-}" && -n "${CXX:-}" ]] || die "Conda compilers \$CC/\$CXX not set after activation."
ok "Toolchain: ${BOLD}$(basename "$CXX")${RST}  ($("$CXX" -dumpversion))"
ok "Cantera + fmt + Boost headers staged in ${DIM}\$CONDA_PREFIX/include${RST}"

# ----------------------------------------------------------------------------
#  3. configure (clean, out-of-source, conda toolchain)
# ----------------------------------------------------------------------------
step "Configuring the build (clean slate)"
rm -rf build
run_spin "Running CMake configure" \
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX"
ok "Build files generated in ${BOLD}build/${RST}"

# ----------------------------------------------------------------------------
#  4. build on $CORES cores  (streamed live — watch the eddies fly)
# ----------------------------------------------------------------------------
step "Compiling odt.x on ${BOLD}${CORES}${RST} cores"
echo "${DIM}  (live compiler output below)${RST}"
echo
cmake --build build -j"${CORES}"
echo
ok "Compilation finished."

# ----------------------------------------------------------------------------
#  5. place the binary + sanity check
# ----------------------------------------------------------------------------
step "Finalizing"
BIN=""
for cand in build/src/odt.x build/odt.x; do
  [[ -x "$cand" ]] && BIN="$cand" && break
done
[[ -n "$BIN" ]] || die "Could not find the built odt.x — check the build output above."
ok "Built executable: ${BOLD}${BIN}${RST}"

if [[ -d run ]]; then
  cp -f "$BIN" run/
  ok "Copied to ${BOLD}run/odt.x${RST} (where the run scripts expect it)."
fi

# Confirm the shared libs resolve into the env rather than dangling.
if ldd "$BIN" | grep -Eqi 'cantera|yaml' ; then
  if ldd "$BIN" | grep -Ei 'cantera|yaml' | grep -qi 'not found'; then
    warn "Some libraries did not resolve — only run inside: conda activate ${ENV_NAME}"
  else
    ok "Library check passed (cantera / yaml-cpp resolve into the env)."
  fi
fi

# ----------------------------------------------------------------------------
#  done
# ----------------------------------------------------------------------------
echo
echo "${GRN}${BOLD}  ✔ ODT built successfully.${RST}"
echo
echo "  ${BOLD}To run it:${RST}"
echo "    ${CYN}conda activate ${ENV_NAME}${RST}     ${DIM}# the binary needs the env's libs on its rpath${RST}"
echo "    ${CYN}cd ${ODT_DIR}/run${RST}"
echo "    ${CYN}./odt.x${RST}                  ${DIM}# or one of the run scripts in this folder${RST}"
echo
echo "  ${DIM}Tip: snapshot this environment for reuse elsewhere with${RST}"
echo "    ${DIM}conda env export -n ${ENV_NAME} > odt-env.yml${RST}"
echo
