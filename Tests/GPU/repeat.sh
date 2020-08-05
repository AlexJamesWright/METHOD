#!/bin/sh

cd ../CPU
make test_rk2
./test_rk2
cd ../GPU
make test_rk2
./test_rk2
