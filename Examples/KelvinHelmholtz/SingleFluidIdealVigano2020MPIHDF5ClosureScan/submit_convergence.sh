#!/bin/bash
#SBATCH --job-name=khi-conv
#SBATCH --nodes=1               # iridis requires --nodes; a layout with more ranks than a node has cores (192) needs --nodes=2 on the command line
#SBATCH --ntasks=192            # must equal RX*RY*RZ below (default 8*8*3 for N = 192); override on the command line, see the README
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00         # N = 192: about 1-2 h (estimate); N = 256 on 128 ranks about 6-8 h: submit with --time=12:00:00
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# Convergence run for the 128^3 / 384^3 comparison: the baseline (no closure) of the Vigano et al. KHI at N^3 cells, to t = 14, with
# diagnostics every 20 iterations and primitive frames at the times in times.txt (the 384^3 frame times, see frame_times.py).
# Environment variables (defaults give N = 192 on 192 ranks):
#     N   cells per side          RX RY RZ   process layout (RX*RY*RZ must be --ntasks; each must divide N)          WENO   reconstruction order 5 (default), 7 or 9
#     (directory N<N>_c0 for WENO5, weno<WENO>_N<N>_c0 otherwise; the WENO7 pair for the closure test: N = 256 and 192, section 21 of the design doc)
# Submit from a directory that holds ./main (built from this example) and times.txt, e.g.
#     sbatch                                  submit_convergence.sh                       # N = 192 on 192 ranks (8 8 3)
#     sbatch --ntasks=128 --time=12:00:00 --export=ALL,N=256,RX=8,RY=8,RZ=2 submit_convergence.sh     # N = 256 on 128 ranks
#     sbatch --nodes=1 --ntasks=128 --time=12:00:00 --export=ALL,N=256,RX=8,RY=8,RZ=2,WENO=7 submit_convergence.sh     # WENO7 at 256^3 (estimate 6-7 h)
#     sbatch --nodes=1 --ntasks=192 --time=06:00:00 --export=ALL,WENO=7 submit_convergence.sh                             # WENO7 at 192^3 (8 8 3; estimate 1-2 h)
#     sbatch --nodes=2 --ntasks=256 --time=08:00:00 --export=ALL,N=256,RX=8,RY=8,RZ=4 submit_convergence.sh     # N = 256 on 256 ranks
# --nodes is mandatory on iridis (the submit filter refuses a script without it) and is 1 here; 256 ranks do not fit on one node (an earlier attempt with --nodes=1 and
# --ntasks=256 failed with "Requested partition configuration not available now"): override with --nodes=2 on the command line.
# Build this ./main in its OWN directory: the queued scan tasks use the ./main of the scan directory.
set -euo pipefail

N="${N:-192}"; RX="${RX:-8}"; RY="${RY:-8}"; RZ="${RZ:-3}"; WENO="${WENO:-5}"
[ "$SLURM_NTASKS" -eq $((RX * RY * RZ)) ] || { echo "ERROR: ntasks $SLURM_NTASKS != RX*RY*RZ = $((RX * RY * RZ))"; exit 1; }

source "$HOME/METHOD/modules.sh"
cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
[ -f times.txt ] || { echo "ERROR: times.txt not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

if [ "$WENO" = "5" ]; then RUN="N${N}_c0"; else RUN="weno${WENO}_N${N}_c0"; fi
mkdir -p "$RUN" && cd "$RUN"
echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Run: $RUN ($N^3, layout $RX $RY $RZ)  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ../main 14.0 0 0 0 20 ../times.txt "$N" "$RX" "$RY" "$RZ" "$WENO" > run.log 2>&1
echo "End: $(date)"
