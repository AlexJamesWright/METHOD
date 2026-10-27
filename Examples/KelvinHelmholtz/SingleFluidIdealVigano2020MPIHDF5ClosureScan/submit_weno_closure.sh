#!/bin/bash
#SBATCH --job-name=khi128-wenoc
#SBATCH --nodes=1
#SBATCH --ntasks=64             # must equal RX*RY*RZ below (4*4*4)
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00         # WENO7 baseline: see its sacct line; WENO9 about 5 h (t = 5.3 after about 2 h); generous
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --array=0-12
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

# The closure with a higher-order WENO reconstruction (design doc sections 17-20): does the coarse run, with more small-scale content, need less c to give the exact
# sub-filter exchange? exchange_at_ekin.py on the logs gives C_eff (0.6 for WENO5: section 19); a prediction is that it falls towards the a-priori 0.12-0.3.
# Task: 0 = WENO7 c = 1/8, 1 = WENO7 c = 0.6, 2 = WENO9 c = 1/8, 3 = WENO9 c = 0.6 (the exchange is proportional to c, so two values check the linearity; WENO9 is not robust here, design doc
# section 20: do not submit 2 and 3); 4-6 = WENO7 noise runs, c = 1e-9, 3e-9, 1e-8 (the closure source is physically nil: independent round-off paths of the WENO7 baseline, for its spread);
# 11 = WENO9 and 12 = WENO11 baselines (c = 0), for the early-CFL test below.
# 7-10 = WENO7 at c = 0.2, the exchange-matching value (C_eff 0.14-0.21, design doc section 22), four runs (the 7th digit of c is perturbed for independent round-off paths) to compare with the noise runs.
# Submit from the scan directory (./main built from the CURRENT example, times.txt). The WENO9 tasks (2, 3) are best submitted when the WENO9 baseline has finished:
#     sbatch --array=0,1 submit_weno_closure.sh      # now (WENO7)
#     sbatch --array=4-6 submit_weno_closure.sh      # the WENO7 noise runs
#     sbatch --array=7-10 submit_weno_closure.sh     # WENO7 at the matched c = 0.2
# No frames (only run.log and data_final, saves 3.5 GB per 128^3 run): sbatch --export=ALL,FRAMES=0 ...
# A gentler start (design doc section 23: the WENO9 failures are a transient at t of about 1-3.5): EARLYCFL (default off) and TEARLY, e.g.
#     sbatch --array=11 --export=ALL,EARLYCFL=0.3,TEARLY=4.0 submit_weno_closure.sh      # directory weno9_c0_early0.3
# A different time integrator (default RK3), e.g. the 4th-order SSPRK(5,4) at the same step (5 stages: 1.7 times the cost), design doc section 39:
#     sbatch --array=4 --export=ALL,INTEGRATOR=4,FRAMES=0 submit_weno_closure.sh      # directory weno7_c0.000000001_rk4
# A different CFL number for the whole run (default 0.9): sbatch --array=11 --export=ALL,CFL=0.7 submit_weno_closure.sh      # directory weno9_c0_cfl0.7
# Other grids (N, and the layout RX RY RZ with --ntasks = RX*RY*RZ), e.g. a coarser run with a larger sub-filter share:
#     sbatch --ntasks=8 --array=1,7 --export=ALL,N=64,RX=2,RY=2,RZ=2 submit_weno_closure.sh      # 64^3 (directory weno7_N64_c<c>)
# A separate coefficient for the stress and Poynting terms (1.5 times that of the EMF: design doc section 63): sbatch --array=7 --export=ALL,INTEGRATOR=4,CFL=0.9,QFACTOR=1.5 submit_weno_closure.sh      # directory weno7_c0.2_cfl0.9_rk4_q1.5
# The same run again with frames, in a new directory: sbatch --array=4 --export=ALL,INTEGRATOR=4,CFL=0.9,TAG=_fr submit_weno_closure.sh      # directory weno7_c0.000000001_cfl0.9_rk4_fr
# The continuous Vigano data, optionally Lorentz boosted (design doc sections 69-71; KHInstabilitySingleFluidVigano2020Boosted): BETAX, BETAY, BETAZ (velocity of the new frame, default 0) and
# CONTINUOUS (1 = the continuous data, the default when any of these is set; 0 = the sgn(y) data); when none is set the original Vigano data. Use the nil closure (task 4: c = 1e-9) with RK4:
#     sbatch --array=4 --export=ALL,INTEGRATOR=4,CFL=0.9,CONTINUOUS=1 submit_weno_closure.sh                 # the unboosted continuous control, directory weno7_c0.000000001_cfl0.9_rk4_cont
#     sbatch --array=4 --export=ALL,INTEGRATOR=4,CFL=0.9,BETAZ=0.9 --time=30:00:00 submit_weno_closure.sh     # boosted along z, directory ..._rk4_bz0.9
# The box sees the dynamics dilated by gamma = 1/sqrt(1 - beta^2): ENDTIME (default 14 gamma when boosted, 14 otherwise) and the frame times (the file TIMESFILE, default ../times.txt, times gamma) are scaled, and
# the run needs about gamma times as many steps: raise --time (gamma = 7.1 at beta = 0.99). The times in ../times.txt are those of the 384^3 reference of the ORIGINAL problem: for the continuous problem give
# the end time and a times file of your own (ENDTIME=, TIMESFILE=).
# Output in weno<order>_c<c>/ ; the log has the usual `diag` and `closure` lines.
set -euo pipefail

ORD=(7 7 9 9 7 7 7 7 7 7 7 9 11)
CC=(0.125 0.6 0.125 0.6 0.000000001 0.000000003 0.00000001 0.2 0.2000001 0.2000002 0.2000003 0 0)
W="${ORD[$SLURM_ARRAY_TASK_ID]}"
C="${CC[$SLURM_ARRAY_TASK_ID]}"
N="${N:-128}"; RX="${RX:-4}"; RY="${RY:-4}"; RZ="${RZ:-4}"
[ "$SLURM_NTASKS" -eq $((RX * RY * RZ)) ] || { echo "ERROR: ntasks $SLURM_NTASKS != RX*RY*RZ"; exit 1; }

source "$HOME/METHOD/modules.sh"
cd "$SLURM_SUBMIT_DIR"
[ -x ./main ] || { echo "ERROR: ./main not found in $PWD"; exit 1; }
[ -f times.txt ] || { echo "ERROR: times.txt not found in $PWD"; exit 1; }
export OMP_NUM_THREADS=1

if [ "$N" = "128" ]; then RUN="weno${W}_c${C}"; else RUN="weno${W}_N${N}_c${C}"; fi
EXTRA=""; if [ -n "${EARLYCFL:-}" ]; then EXTRA="$EARLYCFL ${TEARLY:-4.0}"; RUN="${RUN}_early${EARLYCFL}"; fi
if [ -n "${CFL:-}" ]; then [ -n "$EXTRA" ] || EXTRA="0 0"; EXTRA="$EXTRA $CFL"; RUN="${RUN}_cfl${CFL}"; fi          # CFL: the CFL number for the whole run (design doc section 24)
# INTEGRATOR: 2, 3 (default), 4 (SSPRK(5,4)) or 10 (SSPRK(10,4)) (design doc section 39); the directory gets the suffix _rk<INTEGRATOR>
if [ -n "${INTEGRATOR:-}" ]; then
  [ -n "$EXTRA" ] || EXTRA="0 0"
  [ "$(echo $EXTRA | wc -w)" -ge 3 ] || EXTRA="$EXTRA ${CFL:-0.9}"          # the cfl argument precedes the integrator
  EXTRA="$EXTRA $INTEGRATOR"; RUN="${RUN}_rk${INTEGRATOR}"
fi
# QFACTOR: the factor of the quadratic gradient parts of the stress and Poynting deficits (16th argument; the directory gets the suffix _q<QFACTOR>)
if [ -n "${QFACTOR:-}" ]; then
  [ -n "$EXTRA" ] || EXTRA="0 0"
  [ "$(echo $EXTRA | wc -w)" -ge 3 ] || EXTRA="$EXTRA ${CFL:-0.9}"
  [ "$(echo $EXTRA | wc -w)" -ge 4 ] || EXTRA="$EXTRA 3"
  EXTRA="$EXTRA $QFACTOR"; RUN="${RUN}_q${QFACTOR}"
fi
# Continuous / boosted initial data: arguments 17-20 (betaX betaY betaZ continuous) after earlyCfl tEarly cfl integrator quadFactor
BOOSTED=0; GAMMA=1
if [ -n "${BETAX:-}${BETAY:-}${BETAZ:-}${CONTINUOUS:-}" ]; then
  BOOSTED=1
  BX="${BETAX:-0}"; BY="${BETAY:-0}"; BZ="${BETAZ:-0}"
  [ -n "$EXTRA" ] || EXTRA="0 0"
  [ "$(echo $EXTRA | wc -w)" -ge 3 ] || EXTRA="$EXTRA ${CFL:-0.9}"
  [ "$(echo $EXTRA | wc -w)" -ge 4 ] || EXTRA="$EXTRA 3"
  [ "$(echo $EXTRA | wc -w)" -ge 5 ] || EXTRA="$EXTRA 1"
  EXTRA="$EXTRA $BX $BY $BZ ${CONTINUOUS:-1}"
  GAMMA="$(awk -v a="$BX" -v b="$BY" -v c="$BZ" 'BEGIN { print 1 / sqrt(1 - (a * a + b * b + c * c)) }')"
  SUF=""
  for comp in "x:$BX" "y:$BY" "z:$BZ"; do
    v="${comp#*:}"
    if awk -v v="$v" 'BEGIN { exit !(v != 0) }'; then SUF="${SUF}_b${comp%%:*}${v}"; fi
  done
  [ -n "$SUF" ] || SUF="_cont"
  [ "${CONTINUOUS:-1}" = "1" ] || SUF="${SUF}_sgn"
  RUN="${RUN}${SUF}"
fi
RUN="${RUN}${TAG:-}"          # TAG: a suffix for the directory, to repeat a run without overwriting it (e.g. TAG=_fr for the same run with frames)
mkdir -p "$RUN" && cd "$RUN"
ENDTIME="${ENDTIME:-$(awk -v g="$GAMMA" 'BEGIN { print 14.0 * g }')}"
TIMES="${TIMESFILE:-../times.txt}"
if [ "$BOOSTED" = "1" ] && [ "${FRAMES:-1}" != "0" ]; then          # the frame times of the rest frame times gamma (the box sees the dynamics dilated)
  [ -f "$TIMES" ] || { echo "ERROR: times file $TIMES not found"; exit 1; }
  tr ', ' '\n\n' < "$TIMES" | awk -v g="$GAMMA" 'NF { printf "%s%.10g", (n++ ? "," : ""), $1 * g } END { printf "\n" }' > times_scaled.txt
  TIMES="times_scaled.txt"
fi
[ "${FRAMES:-1}" = "0" ] && TIMES="-"          # FRAMES=0: no output frames (the logs only), for runs of which only the log is analysed
echo "Host: $(hostname)  Tasks: $SLURM_NTASKS  Run: $RUN  Frames: $TIMES  End time: $ENDTIME  gamma: $GAMMA  Start: $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core ../main "$ENDTIME" "$C" 0 0 20 "$TIMES" "$N" "$RX" "$RY" "$RZ" "$W" $EXTRA > run.log 2>&1
echo "End: $(date)"
