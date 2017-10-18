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

def getVarFromLine(line, Nx):
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
        return line // (Nx + 1)


def getXIndexFromLine(line, Nx):
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

    Returns
    -------
        index:
            The x-index of the current line's data.
    """
    return (line%(Nx+1) - 1)



# Function declarations over, access data and plot!

if __name__ == '__main__':

    # Get constants first
    with open('../Data/constants.dat', 'r') as f:
        for i, line in enumerate(f):
            if not i==0:
                line=line.split()
                nx = int(line[0])
                ny = int(line[1])
                Nx = int(line[2])
                Ny = int(line[3])
                xmin = float(line[4])
                xmax = float(line[5])
                ymin = float(line[6])
                ymax = float(line[7])
                endTime = float(line[8])
                cfl = float(line[9])
                Ng = int(line[10])
                gamma = float(line[11])
                sigma = float(line[12])
                Ncons = int(line[13])
                Nprims = int(line[14])
                Naux = int(line[15])
                cp = float(line[16])
                dt = float(line[17])
                t = float(line[18])
                dx = float(line[19])
                dy = float(line[20])


    # Now get primitive variables and store the data in array...
    prims = np.zeros([Nprims, Nx, Ny])
    with open('../Data/primitive.dat', 'r') as f:
        for i, line in enumerate(f):
            # Get primitive var labels
            if i==0:
                labels = line.split()[2:]
            # Get primitive var data
            if not i==0 and i%(Nx+1)!=0:
                temp = line.split()
                for j in range(Ny):
                    prims[getVarFromLine(i, Nx)][getXIndexFromLine(i, Nx)][j] = float(temp[j])

    # Clean up labels (remove the commas)
    cleanLabel = []
    for i in range(len(labels)-1):
        cleanLabel.append(labels[i][:-1])
    cleanLabel.append(labels[-1])

    # Plot the heatmaps
    for i in range(Nprims):
        fig = plt.figure()
        plotPrims = prims[i, Ng:-Ng, Ng:-Ng]
        color = cm.plasma
        surf = plt.imshow(plotPrims, cmap=color, interpolation='nearest')
        plt.title(r'Time Evolution for {}: $t = {}$'.format(cleanLabel[i], t))
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.legend()
        plt.show()
