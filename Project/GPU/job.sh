#!/bin/bash

#SBATCH --ntasks-per-node=1     # Tasks per node
#SBATCH --nodes=4                # Number of nodes requested
#SBATCH --partition=gtx1080
#SBATCH --time=00:30:00

module purge
#module load gcc/6.4.0
module load hdf5/1.10.2/gcc/parallel
module load cuda/10.0

module list
nvidia-smi

make clean
make

echo "OMP NUM THREADS"
echo $OMP_NUM_THREADS

#nvprof ./main 10
#mpirun -np 4 ./main 10

mpirun -np 4 nvprof -o methodProfile.%p.nvprof ./main 10



