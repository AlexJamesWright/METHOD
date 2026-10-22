#!/bin/bash
#SBATCH --job-name=khi128-scan
#SBATCH --nodes=1
#SBATCH --ntasks=64             # must equal nxRanks*nyRanks*nzRanks in main.cc (4*4*4)
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00         # my estimate for one run to t = 14 is 1.5 h (extrapolated from the laptop and the 384^3 rate, not measured)
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --array=0-15
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

# Experiment A: the 128^3 closure scan against the 384^3 reference. One array task per configuration; each writes to its own directory
# under the submit directory. Submit from the directory that holds ./main and times.txt (the output times, from frame_times.py):
#     cd /scratch/ih3/METHOD_simulations/vigano2020_128_scan && sbatch submit_scan.sh
# Configurations (c, density term, limiter); task 0 is the baseline without the closure.
set -euo pipefail

CFG=(
  "0 0 0"                  # 0  baseline
  "0.0833333333 0 0"       # 1  c = 1/12
  "0.1666666667 0 0"       # 2  c = 1/6
  "0.3333333333 0 0"       # 3  c = 1/3
  "0.6666666667 0 0"       # 4  c = 2/3
  "1.0 0 0"                # 5  c = 1
  "2.0 0 0"                # 6  c = 2
  "0.3333333333 1 0"       # 7  c = 1/3 with the density term
  "1.0 0 1"                # 8  c = 1 with the limiter
  "2.0 0 1"                # 9  c = 2 with the limiter
  "0.000001 0 0"           # 10 c = 1e-6: the baseline up to rounding, a different round-off path = the noise floor of the comparisons
  "0.000000001 0 0"        # 11 c = 1e-9  } the same, more samples of the noise: each is an independent round-off path of the baseline, to measure how wide the
  "0.000000003 0 0"        # 12 c = 3e-9  } spread of nominally identical runs is (the scan of c shows differences of up to 3 time units at E_mag = 2e-5
  "0.00000001 0 0"         # 13 c = 1e-8  } that are not ordered in c)
  "0.00000003 0 0"         # 14 c = 3e-8  }
  "0.0000001 0 0"          # 15 c = 1e-7  }
)
read -r C DENS LIM <<< "${CFG[$SLURM_ARRAY_TASK_ID]}"

source "$HOME/METHOD/modules.sh"
cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
[ -f times.txt ] || { echo "ERROR: times.txt not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

RUN="c${C}_d${DENS}_l${LIM}"
mkdir -p "$RUN" && cd "$RUN"
echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Run: $RUN  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ../main 14.0 "$C" "$DENS" "$LIM" 20 ../times.txt > run.log 2>&1
echo "End: $(date)"
