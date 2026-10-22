#!/bin/bash
#SBATCH --job-name=khi128-weno
#SBATCH --nodes=1
#SBATCH --ntasks=64             # must equal RX*RY*RZ below (4*4*4)
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00         # estimate 3-4 h per run (2.5-4 times the WENO5 baseline's 1 h, not measured); generous
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --array=0-1
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

# The 128^3 baseline (no closure) with a higher-order WENO reconstruction, to test the estimate of planning/closure_stage2_design.md section 17:
# task 0 = WENO7 (Ng = 5), task 1 = WENO9 (Ng = 6); Lax-Friedrichs flux splitting as in every other run. Same initial data, times.txt and output as the
# WENO5 baseline c0_d0_l0 of the scan. Submit from the scan directory (the one with ./main built from this example and times.txt):
#     cd /scratch/ih3/METHOD_simulations/vigano2020_128_scan && sbatch submit_weno.sh
# Build ./main from the CURRENT version of this example (it takes the order as its 11th argument); the default order 5 reproduces the existing runs.
# Compare with the WENO5 baseline: time_shift.py against the 384^3 log, compare_logs.py, scan_summary.py, compare_scan.py --match ekin.
set -euo pipefail

ORD=(7 9)
W="${ORD[$SLURM_ARRAY_TASK_ID]}"
N=128; RX=4; RY=4; RZ=4
[ "$SLURM_NTASKS" -eq $((RX * RY * RZ)) ] || { echo "ERROR: ntasks $SLURM_NTASKS != RX*RY*RZ"; exit 1; }

source "$HOME/METHOD/modules.sh"
cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
[ -f times.txt ] || { echo "ERROR: times.txt not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

RUN="weno${W}_N${N}_c0"
mkdir -p "$RUN" && cd "$RUN"
echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Run: $RUN  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ../main 14.0 0 0 0 20 ../times.txt "$N" "$RX" "$RY" "$RZ" "$W" > run.log 2>&1
echo "End: $(date)"
