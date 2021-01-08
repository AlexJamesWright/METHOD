#!/bin/bash

cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomSerial || exit
make
# run serial code (the argument is the random seed)
./main 10
cd ../../../Scripts || exit
# copy Data to DataParallel
bash copyData.sh "DataSerial"
python3 text2hdf.py "DataSerial"

cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomSerialHDF5/ || exit
make
# run parallel code (the argument is the random seed)
./main 10
cp data_serial.hdf5 ../../../Scripts
cd ../../../Scripts || exit

# compare converted serial data and native HDF5 output
python3 compareHDF5.py "DataSerial.converted.hdf5" "data_serial.hdf5"




