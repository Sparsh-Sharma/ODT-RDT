#!/bin/bash
# Test-3 Option-A campaign, one case per submission:
#     sbatch --export=CASE=homogeneousStrain2  slrm_test3_array.sh
#     sbatch --export=CASE=homogeneousStrain2A slrm_test3_array.sh
# 128 array tasks x 8 shifts each = 1024 realizations (shifts 0..1023,
# seeds 22..1045 via seed += shift), written to ../../data/$CASE/data/data_NNNNN.
# Paired seeds across the two cases: submit both with the same array range.

#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-127
#SBATCH -J "test3"
##SBATCH --account=YOUR_CARO_PROJECT      # <-- set your CARO project account
##SBATCH --partition=YOUR_PARTITION      # <-- e.g. the standard CPU partition

# module load ...                        # <-- whatever your CARO build used
export LD_LIBRARY_PATH="$HOME/ct/lib:$LD_LIBRARY_PATH"   # conda cantera/yaml-cpp

set -u
: "${CASE:?submit with --export=CASE=<caseName>}"
RUNS_PER_TASK=8

cd "$SLURM_SUBMIT_DIR/.."                # run/ directory (odt.x lives here)

mkdir -p "../data/$CASE/data" "../data/$CASE/input" "../data/$CASE/runtime"
cp -n "../input/$CASE/"* "../data/$CASE/input/" 2>/dev/null || true

for j in $(seq 0 $((RUNS_PER_TASK-1))); do
    shift_no=$(( SLURM_ARRAY_TASK_ID * RUNS_PER_TASK + j ))
    ./odt.x "$CASE" "$shift_no" > "../data/$CASE/runtime/runlog_$shift_no.txt" 2>&1
done
