#!/bin/bash
#SBATCH --job-name=khi384-run
#SBATCH --nodes=1
#SBATCH --ntasks=192            # must equal nxRanks*nyRanks*nzRanks in main.cc (8*8*3)
#SBATCH --cpus-per-task=1
#SBATCH --time=60:00:00         # the queue limit; my estimate for t = 14 is about 30 h (extrapolated, not measured at 384^3)
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# 384^3 Vigano et al. (2020) KHI to t = 14 with energy diagnostics (log lines starting "diag"), primitive frames every 0.5 time
# units and restartable frames every 2 time units. Submit from the run directory that holds ./main:
#     cd /scratch/ih3/METHOD_simulations/vigano2020_384_run1 && sbatch submit_run.sh
# Optional arguments to ./main: endTime diagIters primFrameIters restartFrameIters closureC densityTerm limiter wenoOrder cfl nGrid nxRanks nyRanks nzRanks (see main.cc;
# for WENO7 pass 0 0 0 7 0.9 and give the restart the same order). If the job hits the time limit,
# continue from the latest restart.checkpoint.<t>.hdf5 with ../SingleFluidIdealVigano2020MPIHDF5Restart (nx*ny*nz ranks as here).
set -euo pipefail

source "$HOME/METHOD/modules.sh"

cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ./main 14.0 100 370 1480
echo "End: $(date)"
