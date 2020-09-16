#!/bin/bash

cd KelvinHelmholtz/SingleFluidIdealRandomSerial
make
# run serial code (the argument is the random seed)
./main 10
cd ../..
# copy Data to DataSerial
bash copyData.sh "DataSerial"


cd KelvinHelmholtz/SingleFluidIdealRandomMPI/
make
# run parallel code (the argument is the random seed)
mpirun -np 4 ./main 10
cd ../..
# copy Data to DataParallel
bash copyData.sh "DataParallel"

# compare DataSerial and DataParallel
python3 compare.py "DataSerial" "DataParallel"




