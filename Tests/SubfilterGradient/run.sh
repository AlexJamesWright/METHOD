#!/bin/bash
# Standalone check of SubfilterGradient (no METHOD build): the real class against a numpy reference on analytic fields, for the
# plain closure, with the density term, and with the backscatter limiter; also prints the conservation sums. Needs g++ and numpy.
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../../Project/CPU"
BUILD="$(mktemp -d)"
cd "$BUILD"
cp "$SRC/Include/subfilterGradient.h" "$SRC/Include/modelExtension.h" "$SRC/Src/subfilterGradient.cc" "$HERE/mock/simData.h" "$HERE/test.cc" "$HERE/reference.py" .
g++ -std=c++11 -O2 -I. test.cc subfilterGradient.cc -o test
for args in "0 0 3.0" "1 0 3.0" "0 1 3.0"; do
  echo "== densityTerm limiter c = $args"
  ./test $args out.bin
  python3 reference.py $args out.bin
done
