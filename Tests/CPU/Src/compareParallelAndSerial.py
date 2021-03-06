#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 13:44:06 2018

@author: alex

Tests the precision of the serial and parallel version of METHOD to within
some tolerance. To execute these tests, `make test` from both the serial and
parallel test directories, then from the parallel test directory run the
following command:
`py.test -v Src/compareParallelAndSerial.py`
"""

# Tolerance we want precision
TOL = 1e-15


import sys
from glob import glob

if __name__=='__main__':
    fromSpyder = '../'
else:
    fromSpyder = ''

sys.path.append(fromSpyder+'../../Project/CPU/Src')

from interactivePlot import InteractivePlot as Plot

class HidePrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None
    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        sys.stdout = self._original_stdout



class CompareParallelAndSerial(object):

    Parallels = []
    Serials = []
    Appendicies = []
    Ncons = []
    Nprims = []
    Naux = []
    nx = []
    ny = []
    nz = []
    Ng = []
    xbounds = []
    ybounds = []
    zbounds = []

    def __init__(self):
        self.getFiles()
        print("Running tests...")

    def getFiles(self):

        # For each file, determine the appendix and use interactivePlot to
        # gather the data
        for i, serfile in enumerate(glob(fromSpyder+"../TestData/CPU/Conserved/*")):
            ext = serfile.find('.dat')
            app = serfile.find('Conserved/cons') + len('Conserved.cons')
            appendix = serfile[app:ext]
            self.Appendicies.append(appendix)
            print("Fetching {} data...".format(appendix))

            with HidePrints():
                self.Serials.append(Plot(fromSpyder+"../TestData/Serial/", appendix))
                self.Parallels.append(Plot(fromSpyder+"../TestData/CPU/", appendix))

            self.Ncons.append(self.Serials[i].c['Ncons'])
            self.Nprims.append(self.Serials[i].c['Nprims'])
            self.Naux.append(self.Serials[i].c['Naux'])
            self.nx.append(self.Serials[i].c['nx'])
            self.ny.append(self.Serials[i].c['ny'])
            self.nz.append(self.Serials[i].c['nz'])
            self.Ng.append(self.Serials[i].c['Ng'])

            # Bounds within arrays which do not include ghost cells
            self.xbounds.append((0, self.nx[-1]))
            self.ybounds.append((0, self.ny[-1]))
            self.zbounds.append((0, self.nz[-1]))

# Instantiate the compare class so we have the data
Compare = CompareParallelAndSerial()


# Test functions

def _compareStateVarArrays(serialArray, parallelArray, Obj, nVars):
   for Nv in range(nVars):
       for i in range(*Compare.xbounds[Obj]):
           for j in range(*Compare.ybounds[Obj]):
               for k in range(*Compare.zbounds[Obj]):
                   try:
                       assert(abs((serialArray[Nv, i, j, k] - parallelArray[Nv, i, j, k]) < TOL))
                   except AssertionError:
                       print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
                       assert(abs((serialArray[Nv, i, j, k] - parallelArray[Nv, i, j, k]) < TOL))


# RK2

# OTVortexPeriodic
def test_ConsEquivalentForRKSplitSrmhdPeriodicOTVSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicOTVSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdPeriodicOTVSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicOTVSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdPeriodicOTVSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicOTVSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

# OTVortexOutflow
def test_ConsEquivalentForRKSplitSrmhdOutflowOTVSF():
  Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowOTVSF')
  Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
  _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdOutflowOTVSF():
  Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowOTVSF')
  Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
  _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdOutflowOTVSF():
  Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowOTVSF')
  Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
  _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])



## KelvinHelmholtzRandomInstabilitySingleFluid

def test_ConsEquivalentForRK2SrmhdOutflowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdOutflowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRK2SrmhdOutflowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdOutflowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRK2SrmhdOutflowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdOutflowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRK2SrmhdPeriodicKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdPeriodicKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRK2SrmhdPeriodicKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdPeriodicKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRK2SrmhdPeriodicKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdPeriodicKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRK2SrmhdFlowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdFlowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRK2SrmhdFlowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdFlowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRK2SrmhdFlowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RK2SrmhdFlowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])


## BrioWuSingleFluid

def test_ConsEquivalentForRK2SrmhdOutflowBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdOutflowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRK2SrmhdOutflowBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdOutflowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRK2SrmhdOutflowBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdOutflowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRK2SrmhdPeriodicBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdPeriodicBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRK2SrmhdPeriodicBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdPeriodicBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRK2SrmhdPeriodicBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdPeriodicBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRK2SrmhdFlowBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdFlowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRK2SrmhdFlowBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdFlowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRK2SrmhdFlowBrioWuSF():
   Obj = Compare.Appendicies.index('RK2SrmhdFlowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

# RKSplit

## KelvinHelmholtzRandomInstabilitySingleFluid

def test_ConsEquivalentForRKSplitSrmhdOutflowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdOutflowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdOutflowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRKSplitSrmhdPeriodicKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdPeriodicKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdPeriodicKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRKSplitSrmhdFlowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdFlowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdFlowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdFlowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdFlowKHRandomInstabilitySF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdFlowKHRandomInstabilitySF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])


## BrioWuSingleFluid

def test_ConsEquivalentForRKSplitSrmhdOutflowBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdOutflowBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdOutflowBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdOutflowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRKSplitSrmhdPeriodicBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdPeriodicBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdPeriodicBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdPeriodicBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])

def test_ConsEquivalentForRKSplitSrmhdFlowBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdFlowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.cons, Parallel.cons, Obj, Compare.Ncons[Obj])

def test_PrimsEquivalentForRKSplitSrmhdFlowBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdFlowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.prims, Parallel.prims, Obj, Compare.Nprims[Obj])

def test_AuxEquivalentForRKSplitSrmhdFlowBrioWuSF():
   Obj = Compare.Appendicies.index('RKSplitSrmhdFlowBrioWuSF')
   Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
   _compareStateVarArrays(Serial.aux, Parallel.aux, Obj, Compare.Naux[Obj])
