#!/bin/bash

cd KelvinHelmholtz/SingleFluidIdealRandomMPI
make
# run serial code (the argument is the random seed)
mpirun -np 4 ./main 10
cd ../..
# copy Data to DataParallel
bash copyData.sh "DataParallel"


cd KelvinHelmholtz/SingleFluidIdealRandomMPIHDF5/
make
# run parallel code (the argument is the random seed)
mpirun -np 4 ./main 10
cd ../..
# copy Data to DataParallelHDF5
bash copyData.sh "DataParallelHDF5"

# compare DataParallel and DataParallelHDF5
python3 compare.py "DataParallel" "DataParallelHDF5"




