#!/bin/bash
#SBATCH --job-name=khi128-wenoc
#SBATCH --nodes=1
#SBATCH --ntasks=64             # must equal RX*RY*RZ below (4*4*4)
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00         # WENO7 baseline: see its sacct line; WENO9 about 5 h (t = 5.3 after about 2 h); generous
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --array=0-6
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

# The closure with a higher-order WENO reconstruction (design doc sections 17-20): does the coarse run, with more small-scale content, need less c to give the exact
# sub-filter exchange? exchange_at_ekin.py on the logs gives C_eff (0.6 for WENO5: section 19); a prediction is that it falls towards the a-priori 0.12-0.3.
# Task: 0 = WENO7 c = 1/8, 1 = WENO7 c = 0.6, 2 = WENO9 c = 1/8, 3 = WENO9 c = 0.6 (the exchange is proportional to c, so two values check the linearity; WENO9 is not robust here, design doc
# section 20: do not submit 2 and 3); 4-6 = WENO7 noise runs, c = 1e-9, 3e-9, 1e-8 (the closure source is physically nil: independent round-off paths of the WENO7 baseline, for its spread).
# Submit from the scan directory (./main built from the CURRENT example, times.txt). The WENO9 tasks (2, 3) are best submitted when the WENO9 baseline has finished:
#     sbatch --array=0,1 submit_weno_closure.sh      # now (WENO7)
#     sbatch --array=4-6 submit_weno_closure.sh      # the WENO7 noise runs
# Output in weno<order>_c<c>/ ; the log has the usual `diag` and `closure` lines.
set -euo pipefail

ORD=(7 7 9 9 7 7 7)
CC=(0.125 0.6 0.125 0.6 0.000000001 0.000000003 0.00000001)
W="${ORD[$SLURM_ARRAY_TASK_ID]}"
C="${CC[$SLURM_ARRAY_TASK_ID]}"
N=128; RX=4; RY=4; RZ=4
[ "$SLURM_NTASKS" -eq $((RX * RY * RZ)) ] || { echo "ERROR: ntasks $SLURM_NTASKS != RX*RY*RZ"; exit 1; }

source "$HOME/METHOD/modules.sh"
cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
[ -f times.txt ] || { echo "ERROR: times.txt not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

RUN="weno${W}_c${C}"
mkdir -p "$RUN" && cd "$RUN"
echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Run: $RUN  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ../main 14.0 "$C" 0 0 20 ../times.txt "$N" "$RX" "$RY" "$RZ" "$W" > run.log 2>&1
echo "End: $(date)"
