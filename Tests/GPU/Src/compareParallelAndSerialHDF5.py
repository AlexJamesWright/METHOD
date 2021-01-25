#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 26 Oct 2020 

@author: ania

Tests the precision of the serial plaintext and HDF5 version of METHOD to within
some tolerance. To execute these tests, `make test` from the Tests/CPU directory
"""

import sys
from glob import glob
from pathlib import Path

from compareHDF5 import compare


def test_compareParallelAndSerialHDF5():
   directory1: Path = Path("../TestData/GPUHDF5/")
   directory2: Path = Path("../TestData/MPIGPUHDF5/")

   print("Running tests...")

   # Double check that the previous steps have actually generated the files we expect
   assert(len(list(directory2.glob("*")))>0)
   assert(len(list(directory1.glob("*")))>0)

   # For each file, determine the appendix and use the CompareHDF5 script 
   for serfile in directory2.glob("*"):
       appendix = serfile.stem
       # TODO -- is this still necessary?
       appendix = appendix.strip('aux')
       file1 = directory1 / (appendix + ".hdf5")
       file2 = directory2 / (appendix + ".hdf5")
       print(file1, file2)
       assert(compare(str(file1), str(file2)))
