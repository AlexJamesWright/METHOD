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
        self.final = Plot(self.DataDirectory + '../../Final/', states=False)
        self.Nframes = len(glob(self.DataDirectory + self.variables[0][:-1] + '*'))
        self.t = []
        self.frame = np.zeros([self.final.c['nx'], self.final.c['ny']])



    def _fetchFrame(self, f, v):
        """
        Given an index, update the current frame data for self.frame.

        Parameters
        ----------
        f : int
            Frame number (indexing from 0)
        v : int
            Variable number (index from self.variables)


        """
        var = self.variables[v][:-1]
        print(f"Fetching frame {f+1}/{self.Nframes} of {var}".format(f+1, self.Nframes))

        with open(self.DataDirectory + var + str(f) + '.dat') as file:
            for i, line in enumerate(file):
                if i==0:
                    if v==0:
                        try:
                            self.t.append(float(line[line.index("t = ") + 4:]))
                        except ValueError:
                            # Sometimes theres no header for some reason. Hack a fix.
                            try:
                                self.t.append(self.t[-1] + (self.t[-1]-self.t[-2]))
                            except IndexError:
                                # If it happens in second iterations
                                try:
                                    self.t.append(self.t[-1] + 0.1*self.final.c['dt'])
                                except IndexError:
                                    # If it is the first iteration
                                    self.t.append(0)

                # Get var data
                else:
                    temp = line.split()
                    self.frame[self.final._getXIndexFromLine(i, self.final.c['nx'], self.final.c['ny'])][self.final._getYIndexFromLine(i, self.final.c['nx'], self.final.c['ny'])] = float(temp[0])

        return self.frame

    def _updateImage(self, f, v, im):
        """
        Given an index, update the given image.

        Parameters
        ----------
        f : int
            Frame number (indexing from 0)
        v : int
            Variable number (index from self.variables)
        im : matplotlib.image.AxesImage
            Image to update data with new frame

        Returns
        -------
        im : matplotlib.image.AxesImage
            New frame image

        Notes
        -----
        This is the function to be called iteratively from
        matplotlib.animation.FuncAnimation. The extra variables (v, im) should
        be passed in with the fargs argument.
        """
        self._fetchFrame(f, v)
        im.set_data(self.frame[:, :].T)
        im.set_interpolation('bicubic')
        im.set_clim(self.vmin, self.vmax)
        return im



    def animate(self, filename, v=0, form='.gif', vmin=0.1, vmax=1.5):
        """
        Build and save an animation of a single variable.

        Parameters
        ----------
        v : int
            Variable number (index from self.variables)
        form : str (optional)
            Format to save file. Currently choose from .gif (default) or .mp4
        vmin : float (optional)
            Lower bound for colorbar. Defaults to 0.1
        vmax : float (optional)
            Upper bound for colorbar. Defaults to 1.5

        Notes
        -----
        Save the image in the same directory as animation.py. Normally this
        will be the Src/ directory.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        ext = [self.final.c['xmin'],
               self.final.c['xmax'],
               self.final.c['ymin'],
               self.final.c['ymax']]

        im = ax.imshow(np.zeros((animClass.final.c['nx'], animClass.final.c['ny'])).T,
                       interpolation='bicubic',
                       animated=True,
                       extent=ext,
                       vmin=vmin,
                       vmax=vmax,
                       cmap=cm.CMRmap,
                       origin='lower')

        fig.tight_layout()

        # Create animation
        ani = animation.FuncAnimation(fig,
                                      animClass._updateImage,
                                      animClass.Nframes,
                                      interval=50,
                                      fargs=(v, im))

        # Save the result
        if form == '.gif':
            writer = animation.writers['pillow'](fps=25)
            ani.save(filename+'.gif', writer=writer, dpi=200)
        else:
            writer = animation.writers['ffmpeg'](fps=25)
            ani.save(filename+'.mp4', writer=writer, dpi=200)


if __name__ == '__main__':

    animClass = Anim()
    animClass.animate('METHODAdvert')
