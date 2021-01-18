#!/bin/bash

cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomMPI || exit
make
# run serial code (the argument is the random seed)
mpirun -np 4 ./main 10
cd ../../../Scripts || exit
# copy Data to DataParallel
bash copyData.sh "DataParallel"
python3 text2hdf.py "DataParallel"

cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomMPIHDF5/ || exit
make
# run parallel code (the argument is the random seed)
mpirun -np 4 ./main 10
cp data_parallel.hdf5 ../../../Scripts
cd ../../../Scripts || exit

# compare converted parallel data and native HDF5 output
python3 compareHDF5.py "DataParallel.converted.hdf5" "data_parallel.hdf5"




