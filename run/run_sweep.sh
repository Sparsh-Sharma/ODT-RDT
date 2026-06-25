#!/bin/bash
###############################################################################
# Resolution / adaption sweep to diagnose the Level-1a TKE energy deficit.
#
# Each configuration runs the homogeneousStrain case into its own directory
#   ../data/hs_<label>/data/data_NNNNN/dmp_NNNNN.dat
# so they can be compared side by side.
#
# Primary knob: dxmax (the cell-coarsening cap, as a fraction of domainLength).
# Baseline dxmax=0.02 lets cells merge to 80x the initial spacing (1/4000),
# aggressively coarsening the white-noise field and removing its energy.
# Walking dxmax (and gDens/ngrd0) down should recover kt toward LRR = 4.028.
#
# Run from the run/ directory:   ./run_sweep.sh        (or  ./run_sweep.sh -r  to rebuild first)
# Then analyse:                  python3 analyze_sweep.py
###############################################################################
set -u
echo "start: $(date)"

inputDir="../input/homogeneousStrain"
yamlName=$(basename "$(ls "$inputDir"/*.yaml 2>/dev/null | head -1)")
[ -z "$yamlName" ] && { echo "ERROR: no .yaml found in $inputDir"; exit 1; }
echo "input file: $yamlName"

#       label     ngrd0  gDens  dxmin     dxmax
configs=(
  "baseline   4000   30     0.0002    0.02"
  "dxmax5e3   4000   30     0.0002    0.005"
  "dxmax1e3   4000   60     0.0002    0.001"
  "dxmax5e4   6000   90     0.0001    0.0005"
  "nearunif   8000   120    0.0001    0.0002"
)

rebuild () {
  echo '*** REBUILDING ***'
  ( cd ../build && make -j8 ) || { echo 'FATAL: build error'; exit 1; }
  echo '*** DONE REBUILDING ***'
}
[ "${1:-}" == "-r" ] && rebuild

for cfg in "${configs[@]}"; do
  read -r label ngrd0 gDens dxmin dxmax <<< "$cfg"
  caseName="hs_${label}"
  echo
  echo "===== $caseName : ngrd0=$ngrd0 gDens=$gDens dxmin=$dxmin dxmax=$dxmax ====="

  rm -rf "../data/$caseName"
  mkdir -p "../data/$caseName/data" "../data/$caseName/input" "../data/$caseName/runtime"
  cp "$inputDir/"* "../data/$caseName/input/" 2>/dev/null

  f="../data/$caseName/input/$yamlName"
  sed -i -E "s|^( *ngrd0:[[:space:]]*).*|\1$ngrd0|" "$f"
  sed -i -E "s|^( *gDens:[[:space:]]*).*|\1$gDens|" "$f"
  sed -i -E "s|^( *dxmin:[[:space:]]*).*|\1$dxmin|" "$f"
  sed -i -E "s|^( *dxmax:[[:space:]]*).*|\1$dxmax|" "$f"

  ./odt.x "$caseName" 0
done

echo
echo "===== all runs done; analyse with:  python3 analyze_sweep.py ====="
echo "end: $(date)"
exit 0
