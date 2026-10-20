#!/bin/bash
# Synthetic 3D tests of the checkpoint-restart readers (serial and MPI), with mock Data/ParallelEnv headers so that only the
# real Project/CPU/Src/{initFuncFromCheckpoint,parallelInitFuncFromCheckpoint}.cc are compiled; no METHOD build is needed.
# Needs MPI-enabled HDF5 (h5pcc) and mpicxx. Run from anywhere:  Tests/CheckpointRestart3D/run.sh
# Built with AddressSanitizer: the pre-fix readers overflow their buffer in 3D and abort; the fixed ones pass.
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../../Project/CPU"
BUILD="$(mktemp -d)"
cd "$BUILD"
# Quote-includes search the including file's directory first, so the readers are copied next to the mock headers
cp "$SRC/Src/parallelInitFuncFromCheckpoint.cc" "$SRC/Include/parallelInitFuncFromCheckpoint.h" \
   "$SRC/Src/initFuncFromCheckpoint.cc" "$SRC/Include/initFuncFromCheckpoint.h" "$HERE"/mock/*.h "$HERE"/test_*.cc .
FLAGS="-std=c++17 -g -O0 -fsanitize=address -I. -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lhdf5_hl -lhdf5"
mpicxx $FLAGS test_serial.cc initFuncFromCheckpoint.cc -o test_serial
mpicxx $FLAGS test_parallel.cc parallelInitFuncFromCheckpoint.cc -o test_parallel
./test_serial
for g in "1 1 1" "2 1 1" "1 3 1" "1 1 2" "2 3 2"; do
  n=$(echo $g | awk '{print $1*$2*$3}')
  mpirun --oversubscribe -n "$n" ./test_parallel $g
done
