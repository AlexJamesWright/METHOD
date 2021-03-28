#!/bin/bash

# This script submits a Southampton Iridis5 batch job for the cpu tests
# in Tests/CPU

#SBATCH --ntasks-per-node=4     # Tasks per node
#SBATCH --nodes=1                # Number of nodes requested
#SBATCH --time=00:10:00          # walltime

module purge
module load gcc/6.4.0
module load python/3.6.4
module load hdf5/1.10.2/gcc/parallel
#module load hdf5/1.10.2/gcc/serial

module list
source ../../venv/bin/activate

# -------------- PARAMETERS USERS NEED TO EDIT -------------------

# Enter absolute path to METHOD/Scripts directory here
SCRIPT_DIR=/absolute/path/to/method/root/Scripts

# -----------------------------------------------------------------

# Let python find the scripts for comparing hdf5 files
export PYTHONPATH=$PYTHONPATH:$SCRIPT_DIR

gcc --version
make clean
make test
