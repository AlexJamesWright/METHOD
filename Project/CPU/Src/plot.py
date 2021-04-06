import h5py 
import glob
import os
import numpy as np
# import yt
from abc import ABC, abstractmethod

class DataBase(ABC):

    def __init__(self, datadir):
        self.datadir = datadir 
        self._getExtension()
        
    def _getExtension(self):
        """
        What format is are the files saved in? If both are present, will default to h5.
        """
        exts = ['.dat', '.hdf5']
        for ext in exts:
            if len(glob.glob(self.datadir + '*' + ext)) > 0:
                self._ext = ext
                return self._ext
        raise RuntimeError(f'Not suitable file formats found in {self.datadir}')

    def _getAllMatchingFiles(self, pattern=''):
        """
        Return a list of all files in datadir with the matching extension or pattern.
        """
        if pattern: self._allFiles = glob.glob(self.datadir + '*' + pattern + '*')
        else:       self._allFiles = glob.glob(self.datadir + '*' + self._ext)
        return self._allFiles


class TimeSeriesData(DataBase):

    def __init__(self, datadir):
        super().__init__(datadir)
        self.__ntimes = None
        self.__vars   = None
        
    @property
    def ntimes(self):
        return self.__ntimes
    
    @property 
    def ts(self):
        return self.__ts
    
    @property 
    def vars(self):
        return self.__vars

    def _getFiles(self, pattern=''):
        """
        This will not actually load the data as there may be many checkpoint files (too many to fit in memory). 
        Instead, we will return data in a lazy way. For now, need to get all the file names that contain 
        'checkpoint' in time order.
        """
        self._getAllMatchingFiles(pattern)
        self._files = sorted([file for file in self._allFiles if 'checkpoint' in file], key=os.path.getmtime)
        self.__ntimes = len(self._files)
        self.__ts = self._getFrameTimes()
        self.__vars = self._getVars()
        return self._files

    def _getFrameTimes(self):
        ts = []
        for file in self._files:
            ckpnt = 'checkpoint.'
            startidx = file.index(ckpnt) + len(ckpnt) 
            endidx   = file.index(self._ext)
            ts.append(float(file[startidx:endidx]))
        return np.array(ts)
    
    def _getVars(self):
        return list(self._loadFileFromIndex(0)['UserDef'].keys())


    def _loadFileFromIndex(self, idx):
        """
        Load a h5py file for a specific integar index.
        """
        return h5py.File(self._files[idx], 'r')

    def _loadFileFromTime(self, t):
        """
        Load a h5py file for a specific time. If a frame doesnt exist for the 
        exact time requested, it will return the file for the frame that is
        closest.
        """
        idx = np.argmin(np.abs(self.ts-t))
        return self._loadFileFromIndex(idx)
        
    def __getitem__(self, pos):
        """
        Return the frames 
        Ordering of indices: variable name, time index, x, y, z
        
        Example
        -------
        
        Return 2D rho data for every time frame
        
        >>> dat['rho'].shape
            (10, 256, 128)
        Return the 5=0 pressure data for the final 5 frames
        
        >>> dat['p', -5:, :, 64].shape
            (5, 256)
        """
        if isinstance(pos, str): pos = (pos, slice(None))
        if isinstance(pos[1], slice):
            s =  pos[1]
            stop  = s.stop  or self.ntimes
            start = s.start or 0
            step  = s.step  or 1
            s = slice(start, stop, step)
            return np.stack([np.array(self._loadFileFromIndex(d)['UserDef'][pos[0]])[pos[2:]] for d in range(s.start, s.stop)], axis=0)
        else:
            return np.array(self._loadFileFromIndex(pos[1])['UserDef'][pos[0]])[pos[2:]]
        

# class SingleTimeData(DataBase):
#     def loadData(self, pattern=''):
#         """
#         Need to load the data that does not contain 'checkpoint' in the file name.
#         """
#         self._getAllMatchingFiles(pattern)
#         self._files = [file for file in self._allFiles if 'checkpoint' not in file]
#         if len(self._files) > 1: raise RuntimeError('Found more than one matching file, try a more specific pattern')
#         return self._files


# class DataContainer:
#     pass


# class PlotMETHOD:
#     pass




if __name__=='__main__':
    d = TimeSeriesData('../')
    d._getFiles()
