#!/bin/bash

cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomSerial || exit
make
# run serial code (the argument is the random seed)
./main 10
cd ../../../Scripts || exit
# copy Data to DataSerial
bash copyData.sh "DataSerial"


cd ../Examples/KelvinHelmholtz/SingleFluidIdealRandomMPI/ || exit
make
# run parallel code (the argument is the random seed)
mpirun -np 4 ./main 10
cd ../../../Scripts || exit
# copy Data to DataParallel
bash copyData.sh "DataParallel"

# compare DataSerial and DataParallel
python3 compareText.py "DataSerial" "DataParallel"




