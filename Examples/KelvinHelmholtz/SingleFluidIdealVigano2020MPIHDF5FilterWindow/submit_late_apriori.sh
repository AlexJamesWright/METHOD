#!/bin/bash
#SBATCH --job-name=late-apriori
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G               # about 8 GB at peak (3.6 GB frames); asked for 24 GB for safety
#SBATCH --time=03:00:00         # a-priori about 1-2 min per frame (about 7 frames); the budget reads every frame, the bulk of the time
#SBATCH --partition=batch       # set to the partition you use
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# A-priori EMF test and resolved-energy budget of the 384^3 run in the developed phase (t >= 7.5), at r = 3 (onto 128^3):
# planning/closure_stage2_design.md section 15. Submit from the directory where the tables should be written:
#     cd /scratch/ih3/METHOD_simulations/vigano2020_128_scan_v2 && sbatch ~/METHOD/Examples/KelvinHelmholtz/SingleFluidIdealVigano2020MPIHDF5FilterWindow/submit_late_apriori.sh
# Optional environment (sbatch --export=ALL,TMIN=9,RUN=/path ...): RUN = the 384^3 run directory, TMIN = first time analysed.
set -euo pipefail

module load conda                  # as in the Lagrangian-filtering submit scripts (the module name used to create the env)
eval "$(conda shell.bash hook)"
conda activate lagfilt

F="$HOME/METHOD/Examples/KelvinHelmholtz/SingleFluidIdealVigano2020MPIHDF5FilterWindow"
RUN="${RUN:-/scratch/ih3/METHOD_simulations/vigano2020_384_run1}"
TMIN="${TMIN:-7.5}"
cd "$SLURM_SUBMIT_DIR"
echo "Host: $(hostname)  python: $(which python)  Start: $(date)"
python "$F/apriori_emf.py" --block 3 --dirs "$RUN" --tmin "$TMIN" --every 2 --out apriori_b3_late.txt
python "$F/energy_budget.py" --block 3 --dirs "$RUN" --tmin "$TMIN" --out budget_b3_late.txt
echo "End: $(date)"
