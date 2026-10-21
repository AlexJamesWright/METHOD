#!/bin/bash
#SBATCH --job-name=khi-conv
#SBATCH --nodes=1
#SBATCH --ntasks=192            # must equal RX*RY*RZ below (default 8*8*3 for N = 192); override on the command line, see the README
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00         # N = 192: about 1-2 h (estimate); N = 256 on 128 ranks about 6-8 h: submit with --time=12:00:00
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# Convergence run for the 128^3 / 384^3 comparison: the baseline (no closure) of the Vigano et al. KHI at N^3 cells, to t = 14, with
# diagnostics every 20 iterations and primitive frames at the times in times.txt (the 384^3 frame times, see frame_times.py).
# Environment variables (defaults give N = 192 on 192 ranks):
#     N   cells per side          RX RY RZ   process layout (RX*RY*RZ must be --ntasks; each must divide N)
# Submit from a directory that holds ./main (built from this example) and times.txt, e.g.
#     sbatch                                  submit_convergence.sh                       # N = 192 on 192 ranks (8 8 3)
#     sbatch --ntasks=128 --time=12:00:00 --export=ALL,N=256,RX=8,RY=8,RZ=2 submit_convergence.sh     # N = 256 on 128 ranks
#     sbatch --ntasks=256 --time=08:00:00 --export=ALL,N=256,RX=8,RY=8,RZ=4 submit_convergence.sh     # N = 256 on 256 ranks
# Build this ./main in its OWN directory: the queued scan tasks use the ./main of the scan directory.
set -euo pipefail

N="${N:-192}"; RX="${RX:-8}"; RY="${RY:-8}"; RZ="${RZ:-3}"
[ "$SLURM_NTASKS" -eq $((RX * RY * RZ)) ] || { echo "ERROR: ntasks $SLURM_NTASKS != RX*RY*RZ = $((RX * RY * RZ))"; exit 1; }

source "$HOME/METHOD/modules.sh"
cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
[ -f times.txt ] || { echo "ERROR: times.txt not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

RUN="N${N}_c0"
mkdir -p "$RUN" && cd "$RUN"
echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Run: $RUN ($N^3, layout $RX $RY $RZ)  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ../main 14.0 0 0 0 20 ../times.txt "$N" "$RX" "$RY" "$RZ" > run.log 2>&1
echo "End: $(date)"
