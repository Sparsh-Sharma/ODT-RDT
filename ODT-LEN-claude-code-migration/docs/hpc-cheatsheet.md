# HPC cheat sheet — DLR CARO & CARA

Quick reference for running ODT realizations on the two DLR clusters. **Both run SLURM**
and are open to all DLR institutes. Fill in the `<…>` placeholders from your DLR HPC
onboarding (the login hostnames, account/project ID, partitions, and module names are
institute-specific and not publicly documented).

## The two systems (public specs)

| | **CARO** | **CARA** |
|---|---|---|
| Location | Göttingen (GWDG / Univ. Göttingen data centre) | TU Dresden (ZIH / Lehmann centre, LZR) |
| Nodes | ~1,364 | ~2,280 (+10 GPU nodes) |
| CPU | AMD EPYC, dual-socket, **64 cores/node** | AMD EPYC, dual-socket, **32 cores/node** (some 64-core) |
| GPU | — | 10 nodes × 4 **NVIDIA A100** |
| Storage | ~8.4 PB HDD | ~16.5 PB HDD + 0.6 PB SSD |
| Scheduler | SLURM | SLURM |
| Managed by | GWDG | ZIH |

**Choosing:** interchangeable for CPU realization sweeps; prefer **CARA** only if you
need GPUs. Keep one as primary to avoid scattering data across both.

## Connect

```bash
ssh caro      # or: ssh cara      (after adding Host blocks to ~/.ssh/config)
```
`~/.ssh/config` template is in `hpc/ssh_config.sample` / the migration guide (Part C).

## SLURM basics (both clusters)

```bash
sinfo -s                    # partitions / node availability
squeue --me                 # my queued/running jobs
sbatch  slrmJob_array.sh    # submit a job / job array
scancel <jobid>             # cancel
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,NNodes   # accounting after a run
```

## Job-script skeleton for ODT realizations

ODT realizations are **independent** → use a **job array**. Adapt from the code's
`slrmJob_array.sh`. Fill placeholders:

```bash
#!/bin/bash
#SBATCH --job-name=odt-<case>
#SBATCH --account=<your-slurm-account>      # DLR project/account ID
#SBATCH --partition=<partition>             # e.g. standard/compute — cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=<n>               # CARO: up to 64, CARA: up to 32
#SBATCH --time=<hh:mm:ss>
#SBATCH --array=0-<Nrlz-1>                  # one array task per realization batch
#SBATCH --output=slurm-%x-%A_%a.out

module purge
module load <compiler> <mpi> <cmake> <cantera>   # ← fill from `module avail`

cd $SLURM_SUBMIT_DIR
srun ./odt.x <inputDir> <caseName>_$SLURM_ARRAY_TASK_ID
```

Keep `caseName` unique per array task so realizations don't overwrite each other
(warned about explicitly in Stephens & Lignell 2021).

## Build ODT on a cluster (one time)

```bash
module load <compiler> <mpi> <cmake> <cantera>
cd code/odt/build
# edit the CMake config to point at the cluster's Cantera install
cmake ..
make            # produces odt.x in run/
```

## Data hygiene

- Raw output → cluster **scratch/work** filesystem, not `$HOME`, not git.
- Pull back only **post-processed** data + figures:
  `rsync -av caro:~/odt/data/<case>/post/ results/<case>/`.
- Record per results folder: `input.yaml`, code commit hash, `nRlz`, seeds.

## Fill-me-in (from DLR onboarding)

- CARO login hostname: `__________`   · CARA login hostname: `__________`
- DLR jump/gateway host (if any): `__________`
- SLURM account/project ID: `__________`
- Partitions / QOS: CARO `__________` · CARA `__________`
- Module names: compiler `____` · MPI `____` · CMake `____` · Cantera `____`
- Scratch/work path convention: `__________`
