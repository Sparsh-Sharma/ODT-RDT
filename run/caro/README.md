# Test-3 Option-A campaign on CARO (1024 realizations per case)

## One-time build on a CARO login node
    cd <repo>
    # cantera + yaml-cpp via micromamba (no root needed):
    curl -sSL -o micromamba https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    chmod +x micromamba && ./micromamba create -y -p $HOME/ct -c conda-forge libcantera-devel yaml-cpp
    mkdir -p build_caro && cd build_caro
    CONDA_PREFIX=$HOME/ct CANTERA_INCLUDE_PATH=$HOME/ct/include cmake -DCMAKE_PREFIX_PATH=$HOME/ct ..
    make -j8 && cp src/odt.x ../run/

## Submit (from run/caro/; set account/partition in the script first)
    sbatch --export=CASE=homogeneousStrain2  slrm_test3_array.sh
    sbatch --export=CASE=homogeneousStrain2A slrm_test3_array.sh

Each case: 128 tasks x 8 serial runs (~1 min each) = 1024 realizations,
~10 min/task, ~4 GB of dumps per case.

## Analysis (login node, needs numpy/scipy/matplotlib)
    python3 post/acoustics/results_test3/compare_optionA.py \
        data/homogeneousStrain2 data/homogeneousStrain2A

Prints the paired table (A_low, A_high, u2^2/2kt vs strain, 95% bootstrap CIs)
and writes fig_optionA_compare.{png,pdf} + optionA_table.txt next to the script.
