#!/bin/bash

cd KelvinHelmholtz/SingleFluidIdealRandomSerial
make
# run serial code (the argument is the random seed)
./main 10
cd ../..
# copy Data to DataParallel
bash copyData.sh "DataSerial"


cd KelvinHelmholtz/SingleFluidIdealRandomSerialHDF5/
make
# run parallel code (the argument is the random seed)
./main 10
cd ../..
# copy Data to DataSerialHDF5
bash copyData.sh "DataSerialHDF5"

# compare DataSerial and DataSerialHDF5
python3 compare.py "DataSerial" "DataSerialHDF5"




