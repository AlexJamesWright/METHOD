#!/bin/bash

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

export PYTHONPATH=$PYTHONPATH:../../Scripts:/home/amb1u19/METHOD_branches/METHOD_dev_hdf5/Scripts

gcc --version
make clean
make test
