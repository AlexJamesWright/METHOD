"""
    Script gathers all data from the animation data folders and stores
    in the anim class ready to run the movie.
"""


from interactivePlot import InteractivePlot as Plot
import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm
import warnings
from copy import deepcopy
from glob import glob

warnings.filterwarnings('ignore', "No labelled objects found. ")

class HidePrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None
    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        sys.stdout = self._original_stdout




# Change this to the relative path to the data you want to plot
# File names must start with e.g. `primitive`
UsualDirectory = '../Data/TimeSeries/UserDef/'

class Anim(object): 
    
    variables = []
    
    def __init__(self, DataDir=None):
        if DataDir==None:
            self.DataDirectory = UsualDirectory
        else:
            self.DataDirectory = DataDir
        with open(self.DataDirectory + 'info') as f:
            for var in f.readlines():
                self.variables.append(var)
        self.final = Plot(self.DataDirectory + '../../Final/')
        self._gatherData()

    def _gatherData(self):
        # Number of frames
        self.vars = np.zeros([len(self.variables), self.final.c['Nx'], self.final.c['Ny'], self.final.c['Nz']])
        self.Nframes = len(glob(self.DataDirectory + self.variables[0][:-1] + '*'))
    
        self.frame = []
        self.t = []
        # Conserved vector first
        print("Gathering data...\nThis may take some time! Found {} frames...".format(self.Nframes))
        for f in range(self.Nframes):
            print("Fetching frame {}/{}".format(f+1, self.Nframes))
            for v, var in enumerate(self.variables):
                with open(self.DataDirectory + self.variables[v][:-1] + str(f) + '.dat') as file:
                    for i, line in enumerate(file):
                        if i==0:
                            if v==0:
                                self.t.append(float(line[line.index("t = ") + 4:]))
                        # Get var data
                        else:
                            temp = line.split()
                            for k in range(self.final.c['Nz']):
                                self.vars[v][self.final._getXIndexFromLine(i, self.final.c['Nx'], self.final.c['Ny'])][self.final._getYIndexFromLine(i, self.final.c['Nx'], self.final.c['Ny'])][k] = float(temp[k])
            self.frame.append(deepcopy(self.vars))
        print("Ready!")
    



if __name__ == '__main__':
    
#     Get data
    animClass = Anim()
    
    # Animate density
    var = 0
    N = 4
    ##### For some reason this doesnt work when inside a function. ########
    ##### Not a disaster atm so will leave it like this            ########
    
    fig = plt.figure()
    
    # ims is a list of lists, each row is a list of artists to draw in the
    # current frame; here we are just animating one artist, the image, in
    # each frame
    ims = []
    for i in range(animClass.Nframes):
        im = plt.imshow(animClass.frame[i][var][N:-N, N:-N, 0], 
                        interpolation='bicubic', 
                        animated=True,
                        cmap=cm.afmhot)
        ims.append([im])
    
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=10)
    
    plt.colorbar()
    plt.show()
#    ani.save('Output.gif', writer='imagemagick', fps=20)
    
    
    