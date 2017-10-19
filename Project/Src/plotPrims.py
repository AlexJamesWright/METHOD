"""
    Script gathers the state vectors stored in the Data directory, specifically
the primitive state vector, and plots the resulting heatmaps. This script is to
be executed as
    `python plotPrims.py`
from the Project/Src directory. If trying to execute from a different location
the relative paths will not match and python will throw an error trying to open
the data file.
"""


import numpy as np
import pylab as plt
from matplotlib import cm

def getVarFromLine(line, Nx, Ny):
    """
    Given the line number that the iterator is on, and the size of the x-domain,
    returns the index of the primitive variable this data belongs to.

    Parameters
    ----------
        line: int
            The line number the file pointer is pointing to. We want to know which
            primitive variable this line's data corresponds to.
        Nx: int
            The total number (incl ghost cells) of domain cells in the x-direction.
        Ny: int
            The total number (incl ghost cells) of domain cells in the y-direction.

    Returns
    -------
        var:
            The primitive variable index of this line's data.

    Other
    -----
        Function will throw a ValueError if trying to get the primitive index
        of the first (zero'th) line.
    """
    if line == 0:
        raise ValueError('Line zero does not contain any data')
    else:
        return ((line-1)//Ny)//Nx


def getXIndexFromLine(line, Nx, Ny):
    """
    Given the line number that the iterator is on, and the size of the x-domain,
    returns the x-index of this line's data.

    Parameters
    ----------
        line: int
            The line number the file pointer is pointing to. We want to know which
            primitive variable this line's data corresponds to.
        Nx: int
            The total number (incl ghost cells) of domain cells in the x-direction.
        Ny: int
            The total number (incl ghost cells) of domain cells in the y-direction.

    Returns
    -------
        index:
            The x-index of the current line's data.
    """
    return ((line-1)//Ny)%Nx

def getYIndexFromLine(line, Nx, Ny):
    """
    Given the line number that the iterator is on, and the size of the y-domain,
    returns the y-index of this line's data.

    Parameters
    ----------
        line: int
            The line number the file pointer is pointing to. We want to know which
            primitive variable this line's data corresponds to.
        Nx: int
            The total number (incl ghost cells) of domain cells in the x-direction.
        Ny: int
            The total number (incl ghost cells) of domain cells in the y-direction.

    Returns
    -------
        index:
            The y-index of the current line's data.
    """
    return (line-1)%Nx



# Function declarations over, access data and plot!

if __name__ == '__main__':

    # Get constants first
    with open('../Data/constants.dat', 'r') as f:
        for i, line in enumerate(f):
            if not i==0:
                line=line.split()
                nx = int(line[0])
                ny = int(line[1])
                nz = int(line[2])
                Nx = int(line[3])
                Ny = int(line[4])
                Nz = int(line[5])
                xmin = float(line[6])
                xmax = float(line[7])
                ymin = float(line[8])
                ymax = float(line[9])
                zmin = float(line[10])
                zmax = float(line[11])
                endTime = float(line[12])
                cfl = float(line[13])
                Ng = int(line[14])
                gamma = float(line[15])
                sigma = float(line[16])
                Ncons = int(line[17])
                Nprims = int(line[18])
                Naux = int(line[19])
                cp = float(line[20])
                dt = float(line[21])
                t = float(line[22])
                dx = float(line[23])
                dy = float(line[24])
                dz = float(line[25])


    # Now get primitive variables and store the data in array...
    prims = np.zeros([Nprims, Nx, Ny, Nz])
    with open('../Data/primitive.dat', 'r') as f:
        for i, line in enumerate(f):
            # Get primitive var labels
            if i==0:
                labels = line.split()[2:]
            # Get primitive var data
            else:
                temp = line.split()
                for k in range(Nz):
                    prims[getVarFromLine(i, Nx, Ny)][getXIndexFromLine(i, Nx, Ny)][getYIndexFromLine(i, Nx, Ny)] = float(temp[k])

    # Clean up labels (remove the commas)
    cleanLabel = []
    for i in range(len(labels)-1):
        cleanLabel.append(labels[i][:-1])
    cleanLabel.append(labels[-1])

    # Plot the heatmaps
    for i in range(Nprims):
        fig = plt.figure()
        plotPrims = prims[i, Ng:-Ng, Ng:-Ng, Ng]
        color = cm.plasma
        surf = plt.imshow(plotPrims, cmap=color, interpolation='nearest')
        plt.title(r'Time Evolution for {}: $t = {}$'.format(cleanLabel[i], t))
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.legend()
        plt.show()
