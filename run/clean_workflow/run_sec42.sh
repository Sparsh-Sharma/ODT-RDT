#!/bin/bash
###############################################################################
# Section 4.2 driver: runs the three canonical cases from one base input,
# stamping the per-run flags so you never hand-edit yaml. Run from run/.
#
#   sp_A1_rdt : eddies OFF, dilatation OFF, inviscid  -> moments baseline
#   sp_A2_dil : eddies OFF, dilatation ON,  inviscid  -> linear-RDT compression
#   sp_B_full : eddies ON,  dilatation ON,  kvisc0=KV -> emergent cascade
#
# Requires the boundary-clamp fix in domain.cc (domainPositionToIndex) for B.
# Base input expected at  ../input/strainSpectrum/input.yaml
#
# Usage:   ./run_sec42.sh           (kvisc0 for B defaults to 1e-5)
#          ./run_sec42.sh 1e-6      (set B's kvisc0)
#          ./run_sec42.sh -r 1e-6   (rebuild first)
###############################################################################
set -u
echo "start: $(date)"

[ "${1:-}" == "-r" ] && { REBUILD=1; shift; } || REBUILD=0
KV="${1:-1.0e-5}"          # kvisc0 for RUN B

baseInput="../input/strainSpectrum"
yamlName=$(basename "$(ls "$baseInput"/*.yaml 2>/dev/null | head -1)")
[ -z "$yamlName" ] && { echo "ERROR: no yaml in $baseInput"; exit 1; }

if [ "$REBUILD" == "1" ]; then
  echo '*** REBUILDING ***'; ( cd ../build && make -j8 ) || { echo 'BUILD FAIL'; exit 1; }
fi

# args: caseName  LnoEddies  Ldilatation  kvisc0
do_run () {
  local caseName="$1" noedd="$2" dil="$3" kv="$4"
  echo; echo "===== $caseName : LnoEddies=$noedd Ldilatation=$dil kvisc0=$kv ====="
  rm -rf "../data/$caseName"
  mkdir -p "../data/$caseName/data" "../data/$caseName/input" "../data/$caseName/runtime"
  cp "$baseInput/"* "../data/$caseName/input/" 2>/dev/null
  local f="../data/$caseName/input/$yamlName"
  sed -i -E "s|^( *LnoEddies:[[:space:]]*).*|\1$noedd|"    "$f"
  sed -i -E "s|^( *Ldilatation:[[:space:]]*).*|\1$dil|"    "$f"
  sed -i -E "s|^( *kvisc0:[[:space:]]*).*|\1$kv|"          "$f"
  local t0=$(date +%s)
  ./odt.x "$caseName" 0 > "../data/$caseName/runtime/stdout.log" 2>&1
  local nd=$(ls "../data/$caseName/data/" 2>/dev/null | wc -l)
  echo "  dumps=$nd  walltime=$(( $(date +%s)-t0 ))s"
  [ "$nd" -lt 5 ] && { echo "  !! crash/early-exit; tail:"; tail -4 "../data/$caseName/runtime/stdout.log" | sed 's/^/     /'; }
}

do_run sp_A1_rdt  true  false 0.0
do_run sp_A2_dil  true  true  0.0
do_run sp_B_full  false true  "$KV"

echo
echo "===== done ====="
echo "Check:"
echo "  A1: python3 analyze_sweep.py      ../data/sp_A1_rdt   (fractions on LRR, deficit ~-0.3%)"
echo "  A2: python3 spectrum_diagnostic.py ../data/sp_A2_dil  (centroid ON exp(-A22 e) -> 7.39)"
echo "  B : python3 spectrum_diagnostic.py ../data/sp_B_full  (centroid ABOVE compression; broadened)"
echo "end: $(date)"
exit 0
