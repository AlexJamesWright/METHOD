"""
This is the scipt used for the magetic reconnection problem analysis. More 
details can be found at arXiv:1906.03150 *****.
"""

from interactivePlot import InteractivePlot as Plot
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os, sys
import seaborn

seaborn.set()

##########################  Function definitions ##############################

class HidePrints(object):
    """
    Neatly and safely hide print statements by redirecting the stdout using the
    "with" keyword.
    """
    def __enter__(self):
        self.orig_out = sys.stdout 
        sys.stdout = open(os.devnull, 'w')
        
    def __exit__(self, _, __, ___):
        sys.stdout.close()
        sys.stdout = self.orig_out
        

def getCurrentDensity(plot):
    """
    Determine the current density for a specified frame of an animClass object.
    
    Parameters
    ----------
    plot : InteractivePlot class
        Containing all the simualtion data
        
    Returns 
    -------
    J : array
        The current density, approximated using Ampere's law.
    """
    idxBx = 5
    idxBy = 6
    N = plot.c['Ng']
    J = (plot.prims[idxBy, 2:, 1:-1, 0] - plot.prims[idxBy, :-2, 1:-1, 0]) / \
    (2*plot.c['dx']) - \
        (plot.prims[idxBx, 1:-1, 2:, 0] - plot.prims[idxBx, 1:-1, :-2, 0]) / \
        (2*plot.c['dy'])
    J = J[N-1:-(N-1), N-1:-(N-1)]
    return J

def gauss(x, a, x0, sigma):
    """
    Definition of the normal Gaussian bell curve.
    
    Parameters
    ----------
    x : array
        x-values of Gaussian curve.
    a : float
        Amplitude.
    x0 : float
        Centre of peak.
    sigma : float
        Standard deviation.
    
    Returns
    -------
    g : array
        Values for the Gaussian bell curve.
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fitGauss(xs, ys):
    """
    Fit a gaussian---defined in the 'gauss' function above---and return the
    parameters that best fit the values given, namely (a, x0, sigma).
    
    Parameters
    ----------
    xs, ys : array
        The x/y values with which to fit a gaussian curve.
        
    Returns 
    -------
    opt : tuple
        The optimum values for (a, x0, sigma) that produce the best fitting 
        Gaussian curve.
    """
    n = len(ys)
    mean = sum(xs*ys) / n
    sigma = sum(ys*(xs-mean)**2)/n
    opt, pcov = curve_fit(gauss, xs, ys, p0=[0, mean, sigma])
    return opt

############################## Import all data ################################
# The data imported here is the final simulation data recorded at time T=100.0
# from the magnetic reconnection problem, with the conductivities given below.
# To run this script successfully, save all data for ideal (SRMHD), resistive
# (SRRMHD) and REGIME in a directory, we will refer to this as 'root'. Data for
# the resistive and REGIME simulations should be split according to the 
# conductivity, and in a directory under that name. E.g. the saved data for one
# of the REGIME simulations should lie in:
#       root/ 
#         \____REGIME/
#                 \_____Sig100/
#                          \____Final/
# where final is the normal directory saved by the 'saveData' class.
    
# Change this as directed above...
root = None

if root is None:
    root = '/media/alex/My Passport/DataToKeep/BigIriData/' + \
            'REGIMETesting/AllData/MagReconnectionPaper'
    
sigs = ['100', '250', '500', '750', '1000']
 
if 'ideal' not in locals():
    print(f"Getting ideal...")
    with HidePrints():
        ideal = Plot(root + '/SRMHD/Final/')
if 'resistives' not in locals():
    resistives = []
    REGIMEs = []
    for sig in sigs:
        print(f"Getting sig = {sig}...")
        with HidePrints():
            resistives.append(Plot(root + '/REGIME/Sig' + sig + '/Final/'))
            REGIMEs.append(Plot(root + '/SRRMHD/Sig' + sig + '/Final/'))
      
        
        
        
###############################################################################
################## Plot the final state of the pressure #######################
###############################################################################
            
if True:       
    fig = plt.figure(figsize=(4*3, 6))
    axs = fig.subplots(2, 3, sharex=True, sharey=True)
    for i, (REGIME, resistive) in enumerate(zip([REGIMEs[0], REGIMEs[1], \
           REGIMEs[-1]], [resistives[0], resistives[1], resistives[-1]])):
        axs[0, i].imshow(resistive.prims[4, 4:-4, 4:-4, 0].T, aspect='auto', 
           interpolation='bicubic', extent=[-12.8, 12.8, -6.4, 6.4], \
           vmin = 0.45, vmax = 1.4)
        cax = axs[1, i].imshow(REGIME.prims[4, 4:-4, 4:-4, 0].T, aspect='auto', 
                 interpolation='bicubic', extent=[-12.8, 12.8, -6.4, 6.4], 
                 vmin = 0.45, vmax = 1.4)
        axs[0, i].annotate(r'$Resistive$', xy=(-10, 5), fontsize=10, \
           color='white')
        axs[1, i].annotate(r'$REGIME$', xy=(-10, 5), fontsize=10, \
           color='white')
        axs[0, i].annotate(r'$\sigma = {}$'.format(int(resistive.c['sigma'])), 
           xy=(-10, -5), fontsize=10, color='white')
        axs[1, i].annotate(r'$\sigma = {}$'.format(int(resistive.c['sigma'])), 
           xy=(-10, -5), fontsize=10, color='white')
    
        axs[0, i].grid(False)
        axs[1, i].grid(False)
        axs[1, i].set_xlabel(r'$x$')
    axs[0, 0].set_ylabel(r'$y$')
    axs[1, 0].set_ylabel(r'$y$')
    fig.tight_layout()
    cbar_ax = fig.add_axes([1, 0.085, 0.03, 0.88])
    cbar = fig.colorbar(cax, cax=cbar_ax)
    cbar.ax.set_ylabel(r'$p$')
    fig.savefig('Figures/FinalStatePressureNoIdealPaper.pdf', \
        extra_artists=(cbar,), bbox_inches='tight', dpi=300)
    fig.show()


###############################################################################
################# Now do the reconnection rate scaling law ####################
###############################################################################
    
if True:
    mid = resistives[0].c['Ny'] // 2
    left = resistives[0].c['Nx'] // 2
    widthsRes = []
    widthsREGIME = []
    eta = []
    frame = -1
    fig = plt.figure(figsize=(12, 4*len(resistives)))
    axs = fig.subplots(len(resistives), 2)
    
    # For each conductivity, find the best fit gaussian of the current density 
    # from ampere's law and plot to ensure is a good fit. Then determine the 
    # reconnection rate from the width of the profile and plot against 
    # conductivity
    
    for i, (REGIME, resistive) in enumerate(zip(REGIMEs, resistives)):
        
        Jres = getCurrentDensity(resistive)
        JREGIME = getCurrentDensity(REGIME)
        xs = resistive.y[4:-4]
        optRes = fitGauss(xs, Jres[left-4])
        optREGIME = fitGauss(xs, JREGIME[left-4])
        widthsRes.append(np.abs(optRes[2]*2*np.sqrt(2*np.log(2))/2))
        widthsREGIME.append(np.abs(optREGIME[2]*2*np.sqrt(2*np.log(2))/2))
        eta.append(1/resistive.c['sigma'])
        
        
        axs[i, 0].plot(resistive.y[4:-4], Jres[left-4, :], 'b-', \
           label='approx')
        axs[i, 0].plot(xs, gauss(xs, *optRes), 'r--',  label='gauss')
        axs[i, 0].set_ylabel(r'$J_z$')
        axs[i, 0].set_xlabel(r'$y$')
        axs[i, 1].plot(resistive.y[4:-4], JREGIME[left-4, :], 'b-', \
           label='approx')
        axs[i, 1].plot(xs, gauss(xs, *optREGIME), 'r--',  label='gauss')
        axs[i, 1].set_ylabel(r'$J_z$')
        axs[i, 1].set_xlabel(r'$y$')
        
        
        jmin = np.min(Jres[left-4, :])
        jmax = np.max(Jres[left-4, :])
        yJ = jmin-0.05*(jmin-jmax)
        axs[i, 0].annotate(r'$Resistive$', xy=(-6, yJ))
        axs[i, 0].annotate(r'$\sigma = {}$'.format(int(resistive.c['sigma'])),\
           xy=(4, yJ))
        jmin = np.min(JREGIME[left-4, :])
        jmax = np.max(JREGIME[left-4, :])
        yJ = jmin-0.05*(jmin-jmax)
        axs[i, 1].annotate(r'$REGIME$', xy=(-6, yJ))
        axs[i, 1].annotate(r'$\sigma = {}$'.format(int(REGIME.c['sigma'])), \
           xy=(4, yJ))
    fig.tight_layout()
#    fig.savefig('Figures/GaussianFitFinalTime.eps', dpi=300)
    fig.show()
        
    
    # Best fit for the reconnection rate
    pRes = np.polyfit(np.log(eta[1:]), np.log(np.abs(widthsRes[1:])), 1)
    pREGIME = np.polyfit(np.log(eta[1:]), np.log(np.abs(widthsREGIME[1:])), 1)
    fig = plt.figure(figsize=(6, 8))
    axs = fig.subplots(2, 1)
    axs[0].loglog(eta, widthsRes, '.')
    axs[0].loglog(eta, np.exp(pRes[1])*eta**pRes[0], \
       label=r'$Gradient \approx {:.2}$'.format(pRes[0]))
    axs[0].set_xlabel(r'$\eta$')
    axs[0].set_ylabel(r'$\lambda / L$')
    axs[0].legend(loc=(0.06, 0.84))
    axs[0].text(2.65*10**-3, 0.082, r'$Resistive$')
    
    axs[1].loglog(eta, widthsREGIME, '.')
    axs[1].loglog(eta, np.exp(pREGIME[1])*eta**pREGIME[0], \
       label=r'$Gradient \approx {:.2}$'.format(pREGIME[0]))
    axs[1].set_xlabel(r'$\eta$')
    axs[1].set_ylabel(r'$\lambda / L$')
    axs[1].legend(loc=(0.06, 0.84))
    axs[1].text(2.7*10**-3, 8.35*10**-2, r'$REGIME$')
    fig.tight_layout()
#    fig.savefig('Figures/ScalingLawFinal.pdf', dpi=300)
    fig.show()



