#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 13:44:06 2018

@author: alex
"""

import sys
from glob import glob

if __name__=='__main__':
    fromSpyder = '../'
else:
    fromSpyder = ''

sys.path.append(fromSpyder+'../../Project/Parallel/Src')

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
    Nx = []
    Ny = []
    Nz = []
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
        for i, serfile in enumerate(glob(fromSpyder+"../TestData/Serial/Conserved/*")):
            ext = serfile.find('.dat')
            app = serfile.find('Conserved/cons') + len('Conserved.cons')
            appendix = serfile[app:ext]
            self.Appendicies.append(appendix)
            print("Fetching {} data...".format(appendix))

            with HidePrints():
                self.Serials.append(Plot(fromSpyder+"../TestData/Serial/", appendix))
                self.Parallels.append(Plot(fromSpyder+"../TestData/Parallel/", appendix))
        
            self.Ncons.append(self.Serials[i].c['Ncons'])
            self.Nprims.append(self.Serials[i].c['Nprims'])
            self.Naux.append(self.Serials[i].c['Naux'])
            self.Nx.append(self.Serials[i].c['Nx'])
            self.Ny.append(self.Serials[i].c['Ny'])
            self.Nz.append(self.Serials[i].c['Nz'])
            self.Ng.append(self.Serials[i].c['Ng'])
        
            self.xbounds.append((self.Ng[-1], self.Nx[-1] - self.Ng[-1]))
            if (self.Ny[-1] > 1):
                self.ybounds.append((self.Ng[-1], self.Ny[-1] - self.Ng[-1]))
            else:
                self.ybounds.append((0, 1))
            if (self.Nz[-1] > 1):
                self.zbounds.append((self.Ng[-1], self.Nz[-1] - self.Ng[-1]))
            else:
                self.zbounds.append((0, 1))



# Instantiate the compare class so we have the data
Compare = CompareParallelAndSerial()
            


# Test functions

# IMEX2
def test_ConsEquivalentForSSP2():
    Obj = Compare.Appendicies.index('SSP2')
    Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
    for Nv in range(Compare.Ncons[Obj]):
        for i in range(*Compare.xbounds[Obj]):
            for j in range(*Compare.ybounds[Obj]):
                for k in range(*Compare.zbounds[Obj]):
                    try:
                        assert(abs((Serial.cons[Nv, i, j, k] - Parallel.cons[Nv, i, j, k]) < 1e-12))
                    except AssertionError:
                        print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
                        assert(abs((Serial.cons[Nv, i, j, k] - Parallel.cons[Nv, i, j, k]) < 1e-12))

def test_PrimsEquivalentForSSP2():
    Obj = Compare.Appendicies.index('SSP2')
    Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
    for Nv in range(Compare.Nprims[Obj]):
        for i in range(*Compare.xbounds[Obj]):
            for j in range(*Compare.ybounds[Obj]):
                for k in range(*Compare.zbounds[Obj]):
                    try:
                        assert(abs((Serial.prims[Nv, i, j, k] - Parallel.prims[Nv, i, j, k]) < 1e-12))
                    except AssertionError:
                        print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
                        assert(abs((Serial.prims[Nv, i, j, k] - Parallel.prims[Nv, i, j, k]) < 1e-12))

def test_AuxEquivalentForSSP2():
    Obj = Compare.Appendicies.index('SSP2')
    Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
    for Nv in range(Compare.Naux[Obj]):
        for i in range(*Compare.xbounds[Obj]):
            for j in range(*Compare.ybounds[Obj]):
                for k in range(*Compare.zbounds[Obj]):
                    try:
                        assert(abs((Serial.aux[Nv, i, j, k] - Parallel.aux[Nv, i, j, k]) < 1e-12))
                    except AssertionError:
                        print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
                        assert(abs((Serial.aux[Nv, i, j, k] - Parallel.aux[Nv, i, j, k]) < 1e-12))

## RK2
#def test_ConsEquivalentForRK2():
#    Obj = Compare.Appendicies.index('RK2')
#    Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
#    for Nv in range(Compare.Ncons[Obj]):
#        for i in range(*Compare.xbounds[Obj]):
#            for j in range(*Compare.ybounds[Obj]):
#                for k in range(*Compare.zbounds[Obj]):
#                    try:
#                        print(Nv, i, j, k)
#                        assert(abs((Serial.cons[Nv, i, j, k] - Parallel.cons[Nv, i, j, k]) < 1e-2))
#                    except AssertionError:
#                        print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
#                        assert(abs((Serial.cons[Nv, i, j, k] - Parallel.cons[Nv, i, j, k]) < 1e-2))
#
#def test_PrimsEquivalentForRK2():
#    Obj = Compare.Appendicies.index('RK2')
#    Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
#    for Nv in range(Compare.Nprims[Obj]):
#        for i in range(*Compare.xbounds[Obj]):
#            for j in range(*Compare.ybounds[Obj]):
#                for k in range(*Compare.zbounds[Obj]):
#                    try:
#                        assert(abs((Serial.prims[Nv, i, j, k] - Parallel.prims[Nv, i, j, k]) < 1e-2))
#                    except AssertionError:
#                        print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
#                        assert(abs((Serial.prims[Nv, i, j, k] - Parallel.prims[Nv, i, j, k]) < 1e-2))
#
#def test_AuxEquivalentForRK2():
#    Obj = Compare.Appendicies.index('RK2')
#    Serial, Parallel = Compare.Serials[Obj], Compare.Parallels[Obj]
#    for Nv in range(Compare.Naux[Obj]):
#        for i in range(*Compare.xbounds[Obj]):
#            for j in range(*Compare.ybounds[Obj]):
#                for k in range(*Compare.zbounds[Obj]):
#                    try:
#                        assert(abs((Serial.aux[Nv, i, j, k] - Parallel.aux[Nv, i, j, k]) < 1e-2))
#                    except AssertionError:
#                        print("Error for (Nv, i, j, k) = ({}, {}, {}, {})".format(Nv, i, j, k))
#                        assert(abs((Serial.aux[Nv, i, j, k] - Parallel.aux[Nv, i, j, k]) < 1e-2))

if __name__=='__main__':
    fromSpyder = '../'
    test_ConsEquivalentForSSP2()
    test_PrimsEquivalentForSSP2()
    test_AuxEquivalentForSSP2()
#    test_ConsEquivalentForRK2()
#    test_PrimsEquivalentForRK2()
#    test_AuxEquivalentForRK2()