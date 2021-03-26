#!/bin/bash

# This script submits a Southampton Iridis5 batch job for the gpu tests
# in Tests/GPU

#SBATCH --ntasks-per-node=2     # Tasks per node
#SBATCH --nodes=1                # Number of nodes requested
#SBATCH --partition=gtx1080
#SBATCH --time=00:10:00

module purge
#module load gcc/6.4.0
module load python/3.6.4
module load hdf5/1.10.2/gcc/parallel
module load cuda/8.0

module list

source ../../venv/bin/activate

# -------------- PARAMETERS USERS NEED TO EDIT -------------------

# Enter absolute path to METHOD/Scripts directory here
SCRIPT_DIR=/absolute/path/to/method/root/Scripts

# -----------------------------------------------------------------

# Let python find the scripts for comparing hdf5 files
export PYTHONPATH=$PYTHONPATH:$SCRIPT_DIR

make clean
make gpu_test

# required for GLIBCXX_3.4.21 module to be available for python
module load gcc/6.4.0

make compare_mpi_test



