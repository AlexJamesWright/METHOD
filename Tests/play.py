#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 10:56:02 2020

@author: alex
"""

import sys
sys.path.append('../Project/GPU/Src')
sys.path.append('../Project/CPU/Src')
from interactivePlotGPU import InteractivePlot as PlotGPU
from interactivePlotCPU import InteractivePlot as PlotCPU


parallel = PlotGPU("TestData/GPU/", "RK2")
#serial   = PlotCPU("TestData/Serial/", "RK2")

pp = parallel.prims
sp = serial.prims


#for sv,  pv in zip(serial.prims, parallel.prims):
#    print(f"{np.sum(np.abs(sv-pv) > 1e-15)}/{30**3} failures")