This directory contains the tests for the GPU-capable code, and as such will only run
if you have an Nvidia graphics card installed.

Before running these tests, make sure you have run the serial tests and have set
the MATCH_SERIAL constant in Project/Parallel/include/timeInt.h to 1.

Once this has been done, run

    make test

and

    py.test -v src/compareParallelAndSerial.py

from this directory
