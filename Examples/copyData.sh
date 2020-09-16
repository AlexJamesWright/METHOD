#!/bin/bash

rm -r $1
cp -r Data $1
rm Data/Final/Auxiliary/*.dat
rm Data/Final/Conserved/*.dat
rm Data/Final/Constants/*.dat
rm Data/Final/Primitive/*.dat
rm Data/Final/Domain/*.dat
rm Data/TimeSeries/Auxiliary/*.dat
rm Data/TimeSeries/Conserved/*.dat
rm Data/TimeSeries/Constants/*.dat
rm Data/TimeSeries/Primitive/*.dat
rm Data/TimeSeries/UserDef/*.dat

