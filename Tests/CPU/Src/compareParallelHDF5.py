#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 26 Oct 2020 

@author: ania

Tests the precision of the parallel plaintext and HDF5 version of METHOD to within
some tolerance. To execute these tests, `make test` from the Tests/CPU directory
"""

import sys
from glob import glob
from pathlib import Path

from compareHDF5 import compare


def test_compareSerialHDF5():
   directory1: Path = Path("../TestData/CPUTextToHDF5/")
   directory2: Path = Path("../TestData/CPUHDF5/")

   print("Running tests...")

   # For each file, determine the appendix and use the CompareHDF5 script 
   for i, serfile in enumerate(directory1.glob("Auxiliary/*")):
       appendix = serfile.stem
       appendix = appendix.strip('aux')
       file1 = directory1 / appendix / ".hdf5"
       file1 = directory2 / appendix / ".hdf5"
       assert(compare(file1, file2))


