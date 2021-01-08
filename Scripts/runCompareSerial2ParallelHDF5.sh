#!/bin/bash

cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomSerialHDF5 || exit
make
# run serial code (the argument is the random seed)
./main 10
cp data_serial.hdf5 ../../../Scripts
cd ../../../Scripts || exit


cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomMPIHDF5/ || exit
make
# run parallel code (the argument is the random seed)
mpirun -np 4 ./main 10
cp data_parallel.hdf5 ../../../Scripts
cd ../../../Scripts || exit

# compare serial and parallel HDF5
python3 compareHDF5.py "data_serial.hdf5" "data_parallel.hdf5"




