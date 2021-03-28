#!/bin/bash

# This script submits a Southampton Iridis5 batch job running on 
# 4 procs for the example in Examples/KelvinHelmholtz/SingleFluidIdealRandomMPIHDF5

# Note the seed after the ./main command would need removing for non-random cases.
# (KelvinHelmholtz/SingleFluidIdealMPI)

#SBATCH --ntasks-per-node=4     # Tasks per node
#SBATCH --nodes=1                # Number of nodes requested
#SBATCH --time=00:10:00          # walltime

module purge
module load gcc/6.4.0
module load hdf5/1.10.2/gcc/parallel
#module load hdf5/1.10.2/gcc/serial

module list

make clean
make

mpirun -np 4 ./main 10

