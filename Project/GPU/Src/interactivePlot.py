"""
    Script gathers the state vectors stored in the Data directory and offers
    functionality to plot various elements.
"""


import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erf
from matplotlib import cm
import warnings
from contextlib import suppress

warnings.filterwarnings('ignore', "No labelled objects found. ")

# Change this to the relative path to the data you want to plot
# File names must start with e.g. `primitive`, anything between this
# and `.dat` should be stored in appendix
# By default, this script will gather data for the final condition of the
# simulation at t=t_end. To gather different data, add arguments to the
# constructor to include the path to the directory and any appendages.
FinalDirectory = '../Data/Final/'
appendix = ''

class InteractivePlot(object):

    def __init__(self, DatDirectory=None, append=None, states=True):
        if DatDirectory is None:
            self.DatDir = FinalDirectory
        else:
            self.DatDir = DatDirectory
        if append is None:
            self.appendix = appendix
        else:
            self.appendix = append
        self.gatherData(states)
        print("Ready!")

    def gatherData(self, states):
        """
        Collects and stores all the data required for plotting the final state of
        the system.
        
        Parameters
        ----------
        states : bool
            Load all of the state arrays. If false, only the constants are
            loaded to save time for animation.

        Notes
        -----
        Stores the following public variables:

            cons : array of float
                (Ncons, nx, ny, nz) Array containing the conserved vector
            consLabels : array of string
                (Ncons,) The labels of the conserved elements
            prims : array of float
                (Nprims, nx, ny, nz) Array containing the primitive vector
            primLabels : array of string
                (Nprims,) The labels of the primitive elements
            aux : array of float
                (Naux, nx, ny, nz) Array containing the auxiliary vector
            auxLabels : array of string
                (Naux,) The labels of the auxiliary elements
            c : dictionary
                Dictionary containing all constant data saved in simData. Access
                elements by typing as an argument the constant you want as a string.
                E.g. to get zmax, enter  -->    c['zmax']
                All links are the same as the constant name in the SimData class.

        """

        # Dictionary to hold constants
        self.c = {}
        c = self.c
        # Get constants first
        print("Fetching constants...")
        with open(self.DatDir + 'Constants/constants' + self.appendix + '.dat', 'r') as f:
            for i, line in enumerate(f):
                if not i==0:
                    line=line.split()
                    c['nx'] = int(line[0])
                    c['ny'] = int(line[1])
                    if c['ny'] == 0:
                        c['ny'] = 1
                    c['nz'] = int(line[2])
                    if c['nz'] == 0:
                        c['nz'] = 1
                    c['Nx'] = int(line[3])
                    c['Ny'] = int(line[4])
                    c['Nz'] = int(line[5])
                    c['xmin'] = float(line[6])
                    c['xmax'] = float(line[7])
                    c['ymin'] = float(line[8])
                    c['ymax'] = float(line[9])
                    c['zmin'] = float(line[10])
                    c['zmax'] = float(line[11])
                    c['endTime'] = float(line[12])
                    c['cfl'] = float(line[13])
                    c['Ng'] = int(line[14])
                    c['gamma'] = float(line[15])
                    c['sigma'] = float(line[16])
                    c['Ncons'] = int(line[17])
                    c['Nprims'] = int(line[18])
                    c['Naux'] = int(line[19])
                    c['cp'] = float(line[20])
                    c['dt'] = float(line[21])
                    c['t'] = float(line[22])
                    c['dx'] = float(line[23])
                    c['dy'] = float(line[24])
                    c['dz'] = float(line[25])

        print("{} conserved vectors".format(c['Ncons']))
        print("{} primitive vectors".format(c['Nprims']))
        print("{} auxiliary vectors".format(c['Naux']))

        if states:
            # Now gather conserved data
            self.cons = np.zeros([c['Ncons'], c['nx'], c['ny'], c['nz']])
            print("Fetching conserved variables...")
            with open(self.DatDir + 'Conserved/cons' + self.appendix + '.dat', 'r') as f:
                for i, line in enumerate(f):
                    # Get cons var labels
                    if i==0:
                        consLabels = line.split()[2:]
                        # Get cons var data
                    else:
                        temp = line.split()
                        for k in range(c['nz']):
                            self.cons[self._getVarFromLine(i, c['nx'], c['ny'])][self._getXIndexFromLine(i, c['nx'], c['ny'])][self._getYIndexFromLine(i, c['nx'], c['ny'])][k] = float(temp[k])
    
    
            # Clean up labels (remove the commas)
            self.cleanConsLabels = []
            for i in range(len(consLabels)-1):
                self.cleanConsLabels.append(consLabels[i][:-1])
            self.cleanConsLabels.append(consLabels[-1])
    
            with suppress(FileNotFoundError):
                # Now get primitive variables if  and store the data in array...
                self.prims = np.zeros([c['Nprims'], c['nx'], c['ny'], c['nz']])
                print("Fetching primitive variables...")
                with open(self.DatDir + 'Primitive/prims' + self.appendix + '.dat', 'r') as f:
                    for i, line in enumerate(f):
                        # Get primitive var labels
                        if i==0:
                            primLabels = line.split()[2:]
                        # Get primitive var data
                        else:
                            temp = line.split()
                            for k in range(c['nz']):
                                self.prims[self._getVarFromLine(i, c['nx'], c['ny'])][self._getXIndexFromLine(i, c['nx'], c['ny'])][self._getYIndexFromLine(i, c['nx'], c['ny'])][k] = float(temp[k])
    
                # Clean up labels (remove the commas)
                self.cleanPrimLabels = []
                for i in range(len(primLabels)-1):
                    self.cleanPrimLabels.append(primLabels[i][:-1])
                self.cleanPrimLabels.append(primLabels[-1])
    
            with suppress(FileNotFoundError):
                # And finally the aux vars if available
                self.aux = np.zeros([c['Naux'], c['nx'], c['ny'], c['nz']])
                print("Fetching auxiliary variables...")
                with open(self.DatDir + 'Auxiliary/aux' + self.appendix +'.dat', 'r') as f:
                    for i, line in enumerate(f):
                        # Get cons var labels
                        if i==0:
                            auxLabels = line.split()[2:]
                        # Get cons var data
                        else:
                            temp = line.split()
                            for k in range(c['nz']):
                                self.aux[self._getVarFromLine(i, c['nx'], c['ny'])][self._getXIndexFromLine(i, c['nx'], c['ny'])][self._getYIndexFromLine(i, c['nx'], c['ny'])][k] = float(temp[k])
    
                # Clean up labels (remove the commas)
                self.cleanAuxLabels = []
                for i in range(len(auxLabels)-1):
                    self.cleanAuxLabels.append(auxLabels[i][:-1])
                self.cleanAuxLabels.append(auxLabels[-1])
            
#            with suppress(FileNotFoundError):
#                # Grab domain data
#                self.x = np.zeros(c['nx'])
#                self.y = np.zeros(c['ny'])
#                self.z = np.zeros(c['nz'])
#                coords = [self.x, self.y, self.z]
#                print("Fetching domain coordinates...")
#                with open(self.DatDir + 'Domain/domain' + self.appendix +'.dat', 'r') as f:
#                    for coord, (i, line) in zip(coords, enumerate(f)):
#                        temp = line.split()
#                        for k, val in enumerate(temp):
#                            coord[k] = float(val)



    def _getVarFromLine(self, line, nx, ny):
        """
        Given the line number that the iterator is on, and the size of the x-domain,
        returns the index of the primitive variable this data belongs to.

        Parameters
        ----------
            line: int
                The line number the file pointer is pointing to. We want to know which
                primitive variable this line's data corresponds to.
            nx: int
                The total number (incl ghost cells) of domain cells in the x-direction.
            ny: int
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
            return ((line-1)//ny)//nx


    def _getXIndexFromLine(self, line, nx, ny):
        """
        Given the line number that the iterator is on, and the size of the x-domain,
        returns the x-index of this line's data.

        Parameters
        ----------
            line: int
                The line number the file pointer is pointing to. We want to know which
                primitive variable this line's data corresponds to.
            nx: int
                The total number (incl ghost cells) of domain cells in the x-direction.
            ny: int
                The total number (incl ghost cells) of domain cells in the y-direction.

        Returns
        -------
            index:
                The x-index of the current line's data.
        """
        return ((line-1)//ny)%nx

    def _getYIndexFromLine(self, line, nx, ny):
        """
        Given the line number that the iterator is on, and the size of the y-domain,
        returns the y-index of this line's data.

        Parameters
        ----------
            line: int
                The line number the file pointer is pointing to. We want to know which
                primitive variable this line's data corresponds to.
            nx: int
                The total number (incl ghost cells)n of domain cells in the x-direction.
            ny: int
                The total number (incl ghost cells) of domain cells in the y-direction.

        Returns
        -------
            index:
                The y-index of the current line's data.
        """
        return (line-1)%ny




    ###############################################################################
    #                             Plotting  Functions                             #
    ###############################################################################




    def plotHeatMaps(self, data='prims', color=None, axis=2):
        """
        Plots the 2D heatmap of the given data. The axes to be plotted can be
        selected via the axis parameter---this corresponds to the axis you want
        to ignore.

        Parameters
        ----------
            data: string
                Describes which variables the user wants to plot. Choose from
                'prims', 'cons', 'aux' or 'primitive', 'conserved' and 'auxiliary'
            color: matplotlib color map
                The colour theme to be plotting in. This can take string arguments
                but best to stick to variants of cm.somecolourscheme
                E.g. cm.magma
            axis: int
                The axis the user wants to ignore.
                (0, 1, 2) = (x, y, z)
        """
        if data=='prims' or data=='primitive':
            data = self.prims
            dataLabels = self.cleanPrimLabels
        elif data=='cons' or data=='conserved':
            data = self.cons
            dataLabels = self.cleanConsLabels
        elif data=='aux' or data=='auxiliary':
            data = self.aux
            data = self.cleanAuxLabels
        else:
            raise ValueError("Variable type not recognised, please try again")
        c = self.c

        for i in range(data.shape[0]):
            fig, ax = plt.subplots(1)
            if (axis == 0):
                plotVars = data[i, c['Nx']//2, :, :]
                axisLabel1 = r'$y$'
                axisLabel2 = r'$z$'
            if (axis == 1):
                plotVars = data[i, :, c['Ny']//2, :]
                axisLabel1 = r'$x$'
                axisLabel2 = r'$z$'
            if (axis == 2):
                plotVars = data[i, :, :, c['Nz']//2]
                axisLabel1 = r'$x$'
                axisLabel2 = r'$y$'

            if color==None:
                color = cm.afmhot
            ext = [self.c['xmin'], self.c['xmax'], self.c['ymin'], self.c['ymax']]
            surf = ax.imshow(plotVars.T, cmap=color, interpolation='bicubic', aspect='auto', origin='lower', extent=ext)
            ax.set_title(r'Time Evolution for {}: $t = {}$'.format(dataLabels[i], c['t']))
            ax.set_xlim([self.c['xmin'], self.c['xmax']])
            ax.set_ylim([self.c['ymin'], self.c['ymax']])
            ax.set_xlabel(axisLabel1)
            ax.set_ylabel(axisLabel2)
            fig.colorbar(surf, shrink=0.5, aspect=5)
            plt.show()
        return ax

    def plotSlice(self, data='prims', axis=0):
        """
        Plots the variation of data in the `axis` direction.

        Parameters
        ----------
            data: string
                Describes which variables the user wants to plot. Choose from
                'prims', 'cons', 'aux' or 'primitive', 'conserved' and 'auxiliary'
            color: matplotlib color map
                The colour theme to be plotting in. This can take string arguments
                but best to stick to variants of cm.somecolourscheme
                E.g. cm.magma
            axis: int, optional
                The axis the user wants to plot in.
                (0, 1, 2) = (x, y, z)
                Defaults to axis=0, x-direction.
        """
        if data=='prims' or data=='primitive':
            data = self.prims
            dataLabels = self.cleanPrimLabels
        elif data=='cons' or data=='conserved':
            data = self.cons
            dataLabels = self.cleanConsLabels
        elif data=='aux' or data=='auxiliary':
            data = self.aux
            dataLabels = self.cleanAuxLabels
        else:
            raise ValueError("Variable type not recognised, please try again")
        c = self.c

        Nx, Ny, Nz = c['Nx'], c['Ny'], c['Nz']

        for i in range(len(data)):
            plt.figure()
            if (axis == 0):
                plotVars = data[i, :, Ny//2, Nz//2]
                axisLabel = r'$x$'
                step = c['dx']
                n = c['nx']
                left, right = c['xmin'], c['xmax']
            if (axis == 1):
                plotVars = data[i, Nx//2, :, Nz//2]
                axisLabel = r'$y$'
                step = c['dy']
                n = c['ny']
                left, right = c['ymin'], c['ymax']
            if (axis == 2):
                plotVars = data[i, Nx//2, Ny//2, :]
                axisLabel = r'$z$'
                step = c['dz']
                n = c['nz']
                left, right = c['zmin'], c['zmax']

            ymin = np.min(plotVars)
            ymax = np.max(plotVars)
            rangeY = ymax - ymin
            ylower = ymin - 0.025 * rangeY
            yupper = ymax + 0.025 * rangeY
            xs = np.linspace(left + step/2, right - step/2, n)
            plt.plot(xs, plotVars, label='{}'.format(dataLabels[i]))
            plt.title(r'Time Evolution for {}: $t = {}$'.format(dataLabels[i], c['t']))
            plt.xlabel(axisLabel)
            plt.ylabel(r'$q_{}(x)$'.format(i+1))
            plt.xlim([c['xmin'], c['xmax']])
#            plt.ylim((ylower, yupper))
            plt.legend(loc='lower center', fontsize=10)
            plt.show()


    def plotTwoFluidSlice(self):
        """
        Plots the variation of total data in the x-direction of the two fluids.

        """

        c = self.c
        Ny, Nz = c['Ny'], c['Nz']

        rho = self.prims[0, :, Ny//2, Nz//2] + self.prims[5, :, Ny//2, Nz//2]
        p   = self.prims[4, :, Ny//2, Nz//2] + self.prims[9, :, Ny//2, Nz//2]
        var = [rho, *self.aux[31:34, :, Ny//2, Nz//2], p, *self.prims[10:, :, Ny//2, Nz//2]]
        varLab = [r'$\rho$', r'$u_x$', r'$u_y$', r'$u_z$', r'$p$', r'$B_x$', r'$B_y$', r'$B_z$', r'$E_x$', r'$E_y$', r'$E_z$']

        xs = np.linspace(c['xmin'] + c['dx']/2, c['xmax'] - c['dx']/2, c['nx'])

        for i, v in enumerate(var):
            plt.figure()
            plt.plot(xs, v)
            plt.title(varLab[i])
            ymin = np.min(v)
            ymax = np.max(v)
            rangeY = ymax - ymin
            ylower = ymin - 0.025 * rangeY
            yupper = ymax + 0.025 * rangeY
            plt.title(r'Time Evolution for {}: $t = {}$'.format(varLab[i], c['t']))
            plt.xlabel(r'$x$')
            plt.ylabel(r'$q_{}(x)$'.format(i+1))
            plt.ylim((ylower, yupper))
            plt.xlim([c['xmin'], c['xmax']])
            plt.legend(loc='lower center', fontsize=10)
            plt.show()

    def plotTwoFluidCurrentSheetAgainstExact(self):
        """
        The current sheet has an analytical solution for the y-direction magnetic
        field. This is plotted against the given B-field.
        """
        By = self.cons[11]
        c = self.c
        plt.figure()
        xs = np.linspace(c['xmin'], c['xmax'], c['nx'])
        exact = np.sign(xs)*erf(0.5 * np.sqrt(c['sigma'] * xs ** 2 / (c['t']+1)))
        plt.plot(xs, By[:, 0, 0], label='Numerical')
        plt.plot(xs, exact, label='Exact')
        plt.xlim([c['xmin'], c['xmax']])
        plt.ylim([-1.2, 1.2])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$B_y$')
        plt.title(r'Comparison of exact and numerical $B_y$ at $t={:.4f}$'.format(c['t']+1))
        plt.legend(loc='upper left')
        plt.show()
        #return np.linalg.norm(exact - By[:, 0, 0])


    def plotSingleFluidCurrentSheetAgainstExact(self, direction=0):
        """
        The current sheet has an analytical solution for the y-direction magnetic
        field. This is plotted against the given B-field.
        """
        c = self.c
        plt.figure()
        nx = self.c['Nx'] // 2
        ny = self.c['Ny'] // 2
        nz = self.c['Nz'] // 2

        if direction == 0:
            B = self.cons[6, :, ny, nz]
            x = np.linspace(c['xmin'], c['xmax'], c['nx'])
        elif direction == 1:
            B = self.cons[7, nx, :, nz]
            x = np.linspace(c['ymin'], c['ymax'], c['ny'])
        else:
            B = self.cons[5, nx, ny, :]
            x = np.linspace(c['zmin'], c['zmax'], c['nz'])

        exact = np.sign(x)*erf(0.5 * np.sqrt(c['sigma'] * x ** 2 / (c['t']+1)))
        initial = np.sign(x)*erf(0.5 * np.sqrt(c['sigma'] * x ** 2 ))
        plt.plot(x, B, label='Numerical')
        plt.plot(x, exact, 'k--', label='Exact')
        plt.plot(x, initial, label='Initial')
        plt.xlim([c['xmin'], c['xmax']])
        plt.ylim([-1.2, 1.2])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$B_y$')
        plt.title(r'Comparison of exact and numerical $B_y$ at $t={:.4f}$'.format(c['t']+1))
        plt.legend(loc='upper left')
        plt.show()

    def plotTwoFluidCPAlfvenWaveAgainstExact(self):
        """
        The cirularly polarized alfven wave has an exact solution, see Amano 2016
        for details. This method plots all non-trivial prims against their exact
        values for case 3.
        """

        rho1, vx1, vy1, vz1, p1, rho2, vx2, vy2, vz2, p2, Bx, By, Bz, Ex, Ey, Ez = self.prims[:]
        c = self.c
        xs = np.linspace(c['xmin'], c['xmax'], c['nx'])
        t = c['t']

        h = 1.04
        B0 = h
        omegaBar1 = -np.sqrt(1.04)
        omegaBar2 = -omegaBar1
        kx = 1.0/4.0

        omega = 5.63803828148e-1
        Wp = 5.19940020571e-6 + 1
        We = 6.68453076522e-5 + 1
        xsi = 0.01

        U1 = -xsi * omega * omegaBar1 / (kx * (omega + omegaBar1 * We))
        U2 = -xsi * omega * omegaBar2 / (kx * (omega + omegaBar2 * Wp))

        phi = kx * xs - omega * t

        BySol = xsi * B0 * np.cos(phi)
        BzSol = -xsi * B0 * np.sin(phi)
        EySol = -(omega/kx)*xsi*B0*np.sin(phi)
        EzSol = -(omega/kx)*xsi*B0*np.cos(phi)
        vy1sol = U1 * np.cos(phi)
        vz1sol = -U1 * np.sin(phi)
        vy2sol = U2 * np.cos(phi)
        vz2sol = -U2 * np.sin(phi)

        # Bx
        BxSol = np.zeros_like(BySol)
        BxSol[:] = B0
        plt.figure()
        plt.plot(xs, Bx[:, 0, 0], label='Numerical')
        plt.plot(xs, BxSol, '--', label='Exact')
        plt.title(r'Exact comparison for $B_x$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # By
        plt.figure()
        plt.plot(xs, By[:, 0, 0], label='Numerical')
        plt.plot(xs, BySol, '--', label='Exact')
        plt.title(r'Exact comparison for $B_y$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # By
        plt.figure()
        plt.plot(xs, Bz[:, 0, 0], label='Numerical')
        plt.plot(xs, BzSol, '--', label='Exact')
        plt.title(r'Exact comparison for $B_z$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # Ex
        plt.figure()
        plt.plot(xs, Ex[:, 0, 0], label='Numerical')
        plt.plot(xs, np.zeros_like(xs), '--', label='Exact')
        plt.title(r'Exact comparison for $E_x$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        minn = min(np.min(Ex), 0)
        maxx = max(np.max(Ex), 0)
        sep = maxx - minn
        plt.ylim([minn-0.1*sep, maxx+0.1*sep])
        plt.legend()
        # Ey
        plt.figure()
        plt.plot(xs, Ey[:, 0, 0], label='Numerical')
        plt.plot(xs, EySol, '--', label='Exact')
        plt.title(r'Exact comparison for $E_y$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # Ez
        plt.figure()
        plt.plot(xs, Ez[:, 0, 0], label='Numerical')
        plt.plot(xs, EzSol, '--', label='Exact')
        plt.title(r'Exact comparison for $E_z$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # vx1
        plt.figure()
        plt.plot(xs, vx1[:, 0, 0], label='Numerical')
        plt.plot(xs, np.zeros_like(xs), '--', label='Exact')
        plt.title(r'Exact comparison for $v_x1$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        minn = min(np.min(vx1), 0)
        maxx = max(np.max(vx1), 0)
        sep = maxx - minn
        plt.ylim([minn-0.1*sep, maxx+0.1*sep])
        plt.legend()
        # vy1
        plt.figure()
        plt.plot(xs, vy1[:, 0, 0], label='Numerical')
        plt.plot(xs, vy1sol, '--', label='Exact')
        plt.title(r'Exact comparison for $v_y1$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # vz1
        plt.figure()
        plt.plot(xs, vz1[:, 0, 0], label='Numerical')
        plt.plot(xs, vz1sol, '--', label='Exact')
        plt.title(r'Exact comparison for $v_z1$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # vx2
        plt.figure()
        plt.plot(xs, vx2[:, 0, 0], label='Numerical')
        plt.plot(xs, np.zeros_like(xs), '--', label='Exact')
        plt.title(r'Exact comparison for $v_x2$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        minn = min(np.min(vx2), 0)
        maxx = max(np.max(vx2), 0)
        sep = maxx - minn
        plt.ylim([minn-0.1*sep, maxx+0.1*sep])
        plt.legend()
        # vy2
        plt.figure()
        plt.plot(xs, vy2[:, 0, 0], label='Numerical')
        plt.plot(xs, vy2sol, '--', label='Exact')
        plt.title(r'Exact comparison for $v_y2$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()
        # vz2
        plt.figure()
        plt.plot(xs, vz2[:, 0, 0], label='Numerical')
        plt.plot(xs, vz2sol, '--', label='Exact')
        plt.title(r'Exact comparison for $v_z2$ at $t={}$'.format(t))
        plt.xlim([c['xmin'], c['xmax']])
        plt.legend()



    def plot2DBrioWu(self, diag=0):
        """
        Plots the main diagonal of the 2D Brio-Wu problem

        Parameters
        ----------
        diag : int
            The diagonal to plot the slice
        """

        nx = self.c['nx']
#        Ny = self.c['Ny']
        midZ = self.c['Nz'] // 2
        Ng = self.c['Ng']

        if diag == 0:
            LB = -Ng
            RB = Ng
            step = -1
        else:
            LB = Ng
            RB = -Ng
            step = 1


        dens = self.prims[0, :, LB:RB:step, midZ].diagonal()
        vx = self.prims[1, :, LB:RB:step, midZ].diagonal()
        vy = self.prims[2, :, LB:RB:step, midZ].diagonal()


        p = self.prims[4, :, LB:RB:step, midZ].diagonal()
        B = self.prims[5, :, LB:RB:step, midZ].diagonal() / np.sqrt(2) + \
            self.prims[6, :, LB:RB:step, midZ].diagonal() / np.sqrt(2)

        # rho
        plt.figure()
        plt.plot(np.linspace(0, 1, nx), dens)
        plt.ylabel(r'$\rho$')
        plt.xlim([0, 1])
        plt.show()
        # vx
        plt.figure()
        plt.plot(np.linspace(0, 1, nx), vx)
        plt.ylabel(r'$vx$')
        plt.xlim([0, 1])
        plt.show()
        # vy
        plt.figure()
        plt.plot(np.linspace(0, 1, nx), vy)
        plt.ylabel(r'$vy$')
        plt.xlim([0, 1])
        plt.show()
        # v rel
        plt.figure()
        plt.plot(np.linspace(0, 1, nx),(vx-vy)/(1-vx*vy))
        plt.ylabel(r'$v (rel)$')
        plt.xlim([0, 1])
        plt.show()
        # v non-rel
        plt.figure()
        plt.plot(np.linspace(0, 1, nx), vx/np.sqrt(2) - vy/np.sqrt(2))
        plt.ylabel(r'$v (non-rel)$')
        plt.xlim([0, 1])
        plt.show()
        # p
        plt.figure()
        plt.plot(np.linspace(0, 1, nx), p)
        plt.ylabel(r'$p$')
        plt.xlim([0, 1])
        plt.show()
        # B
        plt.figure()
        plt.plot(np.linspace(0, 1, nx), B)
        plt.ylabel(r'$B$')
        plt.xlim([0, 1])
        plt.show()

        return B
    
    def plotAdvectionAgainstInitial(self):
        xs = np.linspace(Plot.c['dx']/2, 1-Plot.c['dx']/2, Plot.c['nx'])
        initialRho = np.ones_like(xs)*0.1
        initialRho += 0.4*np.exp(-(10 * (xs - 0.5))**2)
        
        fig, axs = plt.subplots(2)
        fig.set_size_inches(8, 6)
        axs[0].plot(xs, initialRho, 'k-', linewidth=5, alpha=0.3, label='initial')
        axs[0].plot(xs, Plot.prims[0, :, 0, 0], 'b:', label='rho')
        axs[0].set_xlim(xs[0], xs[-1])
        axs[0].set_xlabel(r'$x$')
        axs[0].set_ylabel(r'$\rho$')
        axs[0].legend()
        
        error = np.abs(initialRho-Plot.prims[0, :, 0, 0])
        errorNorm = np.sum(error)/len(error)
        axs[1].semilogy(xs, error, label=rf'Mean = ${errorNorm:.1e}$')
        axs[1].set_xlabel(r"$x$")
        axs[1].set_ylabel('Error')
        axs[1].set_xlim(xs[0], xs[-1])
        axs[1].legend()
        plt.show()
        
        
# Function declarations over, access data and plot!


if __name__ == '__main__':

    Plot = InteractivePlot()

#    Plot.plotHeatMaps()
    
