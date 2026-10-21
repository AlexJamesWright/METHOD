#!/bin/bash
#SBATCH --job-name=khi384-window
#SBATCH --nodes=1
#SBATCH --ntasks=192            # must equal the rank count of the run that wrote the restart file (8*8*3)
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00         # a few snapshots at 384^3 after a short continuation; adjust
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# Write the filtering window at a chosen time by restarting from a restart frame of the 384^3 run. Usage (from a new run directory
# holding ./main built from this example):
#     RESTART=/scratch/.../vigano2020_384_run1/restart.checkpoint.<t>.hdf5 ENDTIME=<window start> NSNAP=30 STRIDE=1 sbatch --export=ALL submit_window.sh
# ENDTIME must be at least 30 iterations (about 0.04 time units at 384^3) after the restart frame's time.
set -euo pipefail

: "${RESTART:?set RESTART to the restart.checkpoint.<t>.hdf5 to continue from}"
: "${ENDTIME:?set ENDTIME to the window start time}"
NSNAP=${NSNAP:-30}
STRIDE=${STRIDE:-1}

source "$HOME/METHOD/modules.sh"

cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Restart: $RESTART  window start $ENDTIME  nSnap $NSNAP stride $STRIDE  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ./main "$RESTART" "$ENDTIME" 0 data_final "$NSNAP" "$STRIDE" 100
echo "End: $(date)"
