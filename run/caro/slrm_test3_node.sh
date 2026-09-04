#!/bin/bash
# Test-3 campaign, node-packed layout for CARO's medium partition, which
# allocates WHOLE 256-core nodes regardless of --cpus-per-task (sacct shows
# AllocCPUS=256 even for 1-CPU tasks).  So: 8 array tasks per case, each
# running 128 serial odt.x realizations CONCURRENTLY on its node.
#     sbatch --export=CASE=<caseName> slrm_test3_node.sh
# Shifts 0..1023 (seeds 22..1045), written to ../../data/$CASE/data/data_NNNNN.
# Same seed pairing across cases as slrm_test3_array.sh.

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=120G
#SBATCH --array=0-7
#SBATCH -J "test3n"

export LD_LIBRARY_PATH="$HOME/ct/lib:$LD_LIBRARY_PATH"   # no-op if ~/ct absent

set -u
: "${CASE:?submit with --export=CASE=<caseName>}"
RUNS_PER_TASK=128

cd "$SLURM_SUBMIT_DIR/.."                # run/ directory (odt.x lives here)

mkdir -p "../data/$CASE/data" "../data/$CASE/input" "../data/$CASE/runtime"
cp -n "../input/$CASE/"* "../data/$CASE/input/" 2>/dev/null || true

base=$(( SLURM_ARRAY_TASK_ID * RUNS_PER_TASK ))
seq "$base" $(( base + RUNS_PER_TASK - 1 )) | \
    xargs -P "$RUNS_PER_TASK" -I{} \
    sh -c 'exec ./odt.x "$1" "$2" > "../data/$1/runtime/runlog_$2.txt" 2>&1' _ "$CASE" {}
