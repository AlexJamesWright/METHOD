"""
    Script gathers all data from the animation data folders and stores
    in the anim class ready to run the movie.
"""


from interactivePlot import InteractivePlot as Plot
import numpy as np
import os, sys
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm
import warnings
from copy import deepcopy

warnings.filterwarnings('ignore', "No labelled objects found. ")

class HidePrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None
    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        sys.stdout = self._original_stdout




# Change this to the relative path to the data you want to plot
# File names must start with e.g. `primitive`, anything between this
# and `.dat` should be stored in appendix
DataDirectory = '../Data/TimeSeries/'

class Anim(object): 
    def __init__(self):
        self._gatherData()

    def _gatherData(self):
        # Number of frames
        self.Nframes = len(os.listdir(DataDirectory + 'Conserved'))
        
        self.frame = []
        # Conserved vector first
        print("Gathering data...\nThis may take some time! Found {} frames...".format(self.Nframes))
        with HidePrints():
            for i, file in enumerate(os.listdir(DataDirectory + 'Conserved')):
                self.frame.append(deepcopy(Plot(DataDirectory, str(i))))
        print("Ready!")
    



if __name__ == '__main__':
    
    # Get data
    animClass = Anim()
    
    # Animate density
    var = 0         
    
    ##### For some reason this doesnt work when inside a function. ########
    ##### Not a disaster atm so will leave it like this            ########
    
    fig = plt.figure()
    
    # ims is a list of lists, each row is a list of artists to draw in the
    # current frame; here we are just animating one artist, the image, in
    # each frame
    ims = []
    for i in range(animClass.Nframes):
        im = plt.imshow(animClass.frame[i].prims[var, 4:-4, 4:-4, 0], interpolation='bicubic', animated=True)
        ims.append([im])
    
    ani = animation.ArtistAnimation(fig, ims, interval=30, blit=True,
                                    repeat_delay=1000)
    
    plt.colorbar()
    plt.show()
    
    
    