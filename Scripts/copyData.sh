#!/bin/bash

rm -r $1
cp -r ../Examples/Data $1
rm ../Examples/Data/Final/Auxiliary/*.dat
rm ../Examples/Data/Final/Conserved/*.dat
rm ../Examples/Data/Final/Constants/*.dat
rm ../Examples/Data/Final/Primitive/*.dat
rm ../Examples/Data/Final/Domain/*.dat
rm ../Examples/Data/TimeSeries/Auxiliary/*.dat
rm ../Examples/Data/TimeSeries/Conserved/*.dat
rm ../Examples/Data/TimeSeries/Constants/*.dat
rm ../Examples/Data/TimeSeries/Primitive/*.dat
rm ../Examples/Data/TimeSeries/UserDef/*.dat

