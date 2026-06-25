#!/bin/bash
###############################################################################
# kvisc0 (Reynolds-number) sweep for the eddies-ON + dilatation-ON run
# (Level 1b, sec:val-spectrum), PARALLEL version.
#
# Runs each kvisc0 as an independent background process into its own
# ../data/hs2_kv<label>/ directory, with at most MAXJOBS in flight at once.
# Each case logs to its own runtime/stdout.log; a summary line prints as each
# finishes (order may vary). eddies ON + dilatation ON are forced per case.
#
# IMPORTANT: the low-viscosity cases are memory- AND time-heavy (ngrd grows,
# tens of millions of eddy trials). MAXJOBS should respect RAM, not just cores
# -- run `free -g` in another terminal during the first batch; if memory is
# climbing toward full, lower MAXJOBS. Do NOT run this until the dilatation/eddy
# boundary fix (domain.cc clamp) is in and a full eddies+dilatation run has been
# confirmed to complete -- otherwise you just crash N cases at once.
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

MAXJOBS=4            # <-- concurrent runs. Keep <= cores AND watch RAM (free -g).

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

run_one () {
  local label="$1" kv="$2"
  local caseName="hs2_${label}"
  rm -rf "../data/$caseName"
  mkdir -p "../data/$caseName/data" "../data/$caseName/input" "../data/$caseName/runtime"
  cp "$inputDir/"* "../data/$caseName/input/" 2>/dev/null
  local f="../data/$caseName/input/$yamlName"
  sed -i -E "s|^( *kvisc0:[[:space:]]*).*|\1$kv|"          "$f"
  sed -i -E "s|^( *LnoEddies:[[:space:]]*).*|\1false|"     "$f"
  sed -i -E "s|^( *Ldilatation:[[:space:]]*).*|\1true|"    "$f"
  local t0=$(date +%s)
  ./odt.x "$caseName" 0 > "../data/$caseName/runtime/stdout.log" 2>&1
  local rc=$?
  local nd=$(ls "../data/$caseName/data/" 2>/dev/null | wc -l)
  echo "[$caseName] kvisc0=$kv  exit=$rc  dumps=$nd  walltime=$(( $(date +%s)-t0 ))s"
  if [ "$nd" -lt 5 ]; then
      echo "   !! few dumps (crash/early-exit) -- tail of its log:"
      tail -4 "../data/$caseName/runtime/stdout.log" | sed 's/^/      /'
  fi
}

###############################################################################
# launch with a concurrency throttle: at most MAXJOBS background jobs at once.
# 'wait -n' (bash >= 4.3, present on WSL Ubuntu) returns when ANY job finishes,
# freeing a slot for the next case.
###############################################################################
for cfg in "${configs[@]}"; do
  read -r label kv <<< "$cfg"
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done
  echo ">>> launching hs2_${label} (kvisc0=${kv})"
  run_one "$label" "$kv" &
done

wait    # block until every case has finished

echo
echo "===== all runs complete; analyse with:  python3 analyze_kvisc_sweep.py ====="
echo "end: $(date)"
exit 0
