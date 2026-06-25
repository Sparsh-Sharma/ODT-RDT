#!/bin/bash
###############################################################################
# kvisc0 (Reynolds-number) sweep for the eddies-ON + dilatation-ON run
# (Level 1b, sec:val-spectrum).
#
# At kvisc0=1e-4 the eddies are so dissipative they destroy the compressed
# small-scale energy faster than the dilatation creates it, so the spectral
# centroid sits BELOW the linear-RDT line exp(-A22 e). Lowering kvisc0 raises
# the Reynolds number, pushes the dissipation scale to higher k, and lets the
# cascade broaden the spectrum. We sweep down by decades to find the window
# where a clean inertial range appears and the centroid climbs toward/above
# the RDT line.
#
# Each kvisc0 runs into ../data/hs2_kv<label>/.  Eddies ON, dilatation ON are
# forced in each copied input regardless of the base file.
#
# Run from run/ :   ./run_kvisc_sweep.sh        (add -r to rebuild first)
# Analyse      :    python3 analyze_kvisc_sweep.py
###############################################################################
set -u
echo "start: $(date)"

inputDir="../input/homogeneousStrain2"
yamlName=$(basename "$(ls "$inputDir"/*.yaml 2>/dev/null | head -1)")
[ -z "$yamlName" ] && { echo "ERROR: no .yaml in $inputDir"; exit 1; }
echo "base input: $inputDir/$yamlName"

#       label    kvisc0
configs=(
  "kv1em4   1.0e-4"
  "kv1em5   1.0e-5"
  "kv1em6   1.0e-6"
  "kv1em7   1.0e-7"
)

rebuild () {
  echo '*** REBUILDING ***'
  ( cd ../build && make -j8 ) || { echo 'FATAL: build error'; exit 1; }
  echo '*** DONE REBUILDING ***'
}
[ "${1:-}" == "-r" ] && rebuild

for cfg in "${configs[@]}"; do
  read -r label kv <<< "$cfg"
  caseName="hs2_${label}"
  echo
  echo "===== $caseName : kvisc0=$kv  (eddies ON, dilatation ON) ====="

  rm -rf "../data/$caseName"
  mkdir -p "../data/$caseName/data" "../data/$caseName/input" "../data/$caseName/runtime"
  cp "$inputDir/"* "../data/$caseName/input/" 2>/dev/null

  f="../data/$caseName/input/$yamlName"
  sed -i -E "s|^( *kvisc0:[[:space:]]*).*|\1$kv|"          "$f"
  sed -i -E "s|^( *LnoEddies:[[:space:]]*).*|\1false|"     "$f"
  sed -i -E "s|^( *Ldilatation:[[:space:]]*).*|\1true|"    "$f"

  t0=$(date +%s)
  ./odt.x "$caseName" 0 > "../data/$caseName/runtime/stdout.log" 2>&1
  rc=$?
  t1=$(date +%s)
  ndump=$(ls "../data/$caseName/data/" 2>/dev/null | wc -l)
  echo "  exit=$rc  walltime=$((t1-t0))s  dumps=$ndump"
  if [ "$ndump" -lt 5 ]; then
      echo "  !! few dumps -- likely crashed/early-exit; tail of log:"
      tail -5 "../data/$caseName/data/$caseName/runtime/stdout.log" 2>/dev/null \
        || tail -5 "../data/$caseName/runtime/stdout.log"
  fi
done

echo
echo "===== done; analyse with:  python3 analyze_kvisc_sweep.py ====="
echo "end: $(date)"
exit 0
