#!/bin/bash
# Standalone serial build and check of KHInstabilitySingleFluidVigano2020Boosted (no HDF5, no MPI). Needs g++ and numpy.
#   ./run.sh build   compile only
#   ./run.sh         compile (if needed) and run the checks
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../../Project/CPU"
BUILD="${BUILD:-$HERE/_build}"
mkdir -p "$BUILD"
cd "$BUILD"
SRCS="simulation.cc srmhd.cc initFunc.cc simData.cc wenoUpwinds.cc weno.cc RK2.cc rkSplit.cc boundaryConds.cc fluxVectorSplitting.cc serialSaveData.cc serialEnv.cc"
CMINPACK="$(ls "$SRC"/CminpackLibrary/Src/*.cc | tr '\n' ' ')"
if [ ! -x test ] || [ "$HERE/test.cc" -nt test ] || [ "$SRC/Src/initFunc.cc" -nt test ] || [ "$SRC/Include/initFunc.h" -nt test ]; then
  OBJS=""
  for s in $SRCS; do
    g++ -std=c++11 -O2 -Wno-unknown-pragmas -I"$SRC/Include" -I"$SRC/CminpackLibrary/Include" -c "$SRC/Src/$s" -o "${s%.cc}.o"
    OBJS="$OBJS ${s%.cc}.o"
  done
  for c in $CMINPACK; do
    g++ -std=c++11 -O2 -I"$SRC/CminpackLibrary/Include" -c "$c" -o "cmin_$(basename "${c%.cc}").o"
    OBJS="$OBJS cmin_$(basename "${c%.cc}").o"
  done
  g++ -std=c++11 -O2 -I"$SRC/Include" -I"$SRC/CminpackLibrary/Include" "$HERE/test.cc" $OBJS -o test
fi
[ "${1:-}" = "build" ] && { echo "built $BUILD/test"; exit 0; }
PY="${PY:-python3}"
for args in "0 0 0 1 1" "0 0 0 0 1" "0 0 0.5 1 1" "0 0 0.9 1 1" "0 0 0.99 1 1" "0.9 0 0 1 1" "0.3 0.4 0.5 1 1" "0 0 0.9 0 1" "0 0 0.9 1 0"; do
  echo "== betaX betaY betaZ continuous mag = $args"
  ./test 32 $args out.bin | sed 's/^/   /'
  $PY "$HERE/reference.py" out.bin $args | sed 's/^/   /'
done
# wider shear layers (the layer width is an argument)
for wargs in "0 0 0 1 1 0.03" "0 0 0.9 1 1 0.03" "0 0 0.8 1 1 0.05"; do
  set -- $wargs
  echo "== betaX betaY betaZ continuous mag layerWidth = $wargs"
  ./test 32 $1 $2 $3 $4 $5 out.bin $6 | sed 's/^/   /'
  $PY "$HERE/reference.py" out.bin $1 $2 $3 $4 $5 $6 | sed 's/^/   /'
done
