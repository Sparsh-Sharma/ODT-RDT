#!/usr/bin/env bash
#
#  rebuild_odt.sh — fast recompile loop for ODT development.
#
#  Does NOT touch conda, the env, or any dependencies. It only:
#    1. Activates the existing 'odt' env (for the toolchain + runtime libs).
#    2. Recompiles changed source, relinks odt.x.
#    3. Copies the fresh binary into run/.
#
#  Two modes:
#    Default (incremental) — recompiles only the files you changed. Fast.
#        bash rebuild_odt.sh
#    Clean — wipes compiled objects + the binary first, full recompile.
#    Use after editing headers widely or CMakeLists.txt.
#        CLEAN=1 bash rebuild_odt.sh
#
#  If you ever change CMakeLists.txt or add/remove source files, prefer CLEAN=1
#  (or re-run the full build_odt.sh) so the build graph is regenerated.
#
set -euo pipefail

ENV_NAME="odt"
CORES=8
MINIFORGE_DIR="${HOME}/miniforge3"
CLEAN="${CLEAN:-0}"

# ----------------------------------------------------------------------------
#  pretty output (same vocabulary as build_odt.sh)
# ----------------------------------------------------------------------------
if [[ -t 1 ]]; then
  BOLD=$'\e[1m'; DIM=$'\e[2m'; RED=$'\e[31m'; GRN=$'\e[32m'
  YLW=$'\e[33m'; BLU=$'\e[34m'; MAG=$'\e[35m'; CYN=$'\e[36m'; RST=$'\e[0m'
else
  BOLD=""; DIM=""; RED=""; GRN=""; YLW=""; BLU=""; MAG=""; CYN=""; RST=""
fi
step()  { echo; echo "${BLU}${BOLD}▸ $*${RST}"; }
info()  { echo "  ${CYN}·${RST} $*"; }
ok()    { echo "  ${GRN}✔${RST} $*"; }
warn()  { echo "  ${YLW}!${RST} $*"; }
die()   { echo; echo "  ${RED}${BOLD}ERROR: $*${RST}" >&2; exit 1; }

trap 'echo; echo "${RED}${BOLD}Rebuild failed.${RST} See the compiler output above."' ERR

START_TS=$(date +%s)
echo
echo "${MAG}${BOLD}  ⟳  ODT rebuild  ${RST}${DIM}(${CLEAN:+clean }$( [[ "$CLEAN" == 1 ]] && echo full || echo incremental ))${RST}"

# ----------------------------------------------------------------------------
#  1. locate the ODT source tree
# ----------------------------------------------------------------------------
step "Locating ODT source"
if [[ -f "CMakeLists.txt" && -d "src" ]]; then
  ODT_DIR="$(pwd)"
elif [[ -f "ODT/CMakeLists.txt" && -d "ODT/src" ]]; then
  ODT_DIR="$(pwd)/ODT"
else
  die "Run this from your ODT tree (the dir with CMakeLists.txt + src/)."
fi
cd "$ODT_DIR"
ok "ODT tree: ${BOLD}${ODT_DIR}${RST}"

# ----------------------------------------------------------------------------
#  2. activate the existing env (do NOT modify it)
# ----------------------------------------------------------------------------
step "Activating '${ENV_NAME}' environment"
CONDA_SH=""
if command -v conda >/dev/null 2>&1; then
  CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
elif [[ -f "${MINIFORGE_DIR}/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="${MINIFORGE_DIR}/etc/profile.d/conda.sh"
else
  die "conda not found. Run build_odt.sh first to set up the environment."
fi
# shellcheck disable=SC1090
source "$CONDA_SH"
conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}" \
  || die "env '${ENV_NAME}' does not exist. Run build_odt.sh first."
conda activate "${ENV_NAME}"
[[ -n "${CC:-}" && -n "${CXX:-}" ]] || die "Toolchain CC/CXX not set after activation."
ok "Toolchain: ${BOLD}$(basename "$CXX")${RST} ($("$CXX" -dumpversion))"

# ----------------------------------------------------------------------------
#  3. ensure the build is configured (only configures if needed)
# ----------------------------------------------------------------------------
if [[ ! -f build/CMakeCache.txt ]]; then
  step "No build/ cache found — configuring once"
  cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
  ok "Configured."
fi

# ----------------------------------------------------------------------------
#  4. clean compiled artifacts if requested (keeps cache/config + deps)
# ----------------------------------------------------------------------------
if [[ "$CLEAN" == "1" ]]; then
  step "Removing compiled objects + binary (deps untouched)"
  cmake --build build --target clean
  rm -f run/odt.x
  ok "Cleaned — next build is a full recompile."
fi

# ----------------------------------------------------------------------------
#  5. compile (incremental unless cleaned) — streamed live
# ----------------------------------------------------------------------------
step "Compiling on ${BOLD}${CORES}${RST} cores"
echo
cmake --build build -j"${CORES}"
echo
ok "Compilation finished."

# ----------------------------------------------------------------------------
#  6. copy fresh binary into run/
# ----------------------------------------------------------------------------
step "Installing fresh odt.x into run/"
BIN=""
for cand in build/src/odt.x build/odt.x; do
  [[ -x "$cand" ]] && BIN="$cand" && break
done
[[ -n "$BIN" ]] || die "odt.x not found after build."
if [[ -d run ]]; then
  cp -f "$BIN" run/
  ok "Copied ${BOLD}${BIN}${RST} → ${BOLD}run/odt.x${RST}"
else
  warn "No run/ directory; binary left at ${BIN}"
fi

# ----------------------------------------------------------------------------
ELAPSED=$(( $(date +%s) - START_TS ))
echo
echo "${GRN}${BOLD}  ✔ Rebuilt in ${ELAPSED}s.${RST}  ${DIM}Run with: conda activate ${ENV_NAME} && cd run && ./odt.x${RST}"
echo
