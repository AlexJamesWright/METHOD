"""
Script used for the power spectrum analysis of the Kelvin-Helmholtz instability
(Figure 6b) in Wright & Hawke (2018).

Resolution for each of model is 512x1024.
"""

from animation import Anim
from interactivePlot import InteractivePlot as Plot
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np

"""
                ###############  NOTE  #################
                #  This script needs to be run twice   #
                #  to execute properly. Dont know why, #
                #  havent tried to figure it out. Go   #
                #  with the flow.                      #
                ########################################
"""

def getFourierTrans(intPlot, u):
    """
    Returns the 1D discrete fourier transform of the variable u along the x-direction
    ready for the power spectrum method.

    Parameters
    ----------
    intPlot : object
        interactivePlot object containing all the simulation data, normally the final instance
    u : ndarray
        Two dimensional array of the variable we want the power spectrum of

    Returns
    -------
    uhat : array (N,)
        Fourier transform of u
    """
    nx, ny = intPlot.c['nx'], intPlot.c['ny']
    NN = nx // 2
    uhat = np.zeros((NN, ny), dtype=np.complex_)

    for k in range(NN):
        for y in range(ny):
            # Sum over all x adding to uhat
            for i in range(nx):
                uhat[k, y] += u[i, y] * np.exp(-(2*np.pi*1j*k*i)/nx)

    return uhat / nx

def getPowerSpectrumSq(intPlot, u):
    """
    Returns the integrated power spectrum of the variable u, up to the Nyquist frequency = nx/2

    Parameters
    ----------
    intPlot : object
        interactivePlot object containing all the simulation data, normally the final instance
    u : ndarray
        Two dimensional array of the variable we want the power spectrum of
    """
    NN = intPlot.c['nx'] // 2
    dy = intPlot.c['dy']
    uhat = getFourierTrans(intPlot, u)
    P = np.zeros(NN)

    for k in range(NN):
        for j in range(intPlot.c['ny']):
            P[k] += (np.absolute(uhat[k, j])**2) * dy

    P = P / np.sum(P)
    return P



def GetKESF(anim, frame):
    """
    Retrieves and computes the kinetic energy density for each frame in a single fluid animation.

    Parameters
    ----------
    anim : object
        animation class containing all user def variables
    frame : Array
        Frame from the animation class containing all user def variables at the time we want

    """
    vx = frame[anim.variables.index("vx\n"), 4:-4, 4:-4, 0]
    vy = frame[anim.variables.index("vy\n"), 4:-4, 4:-4, 0]
    rho = frame[anim.variables.index("rho\n"), 4:-4, 4:-4, 0]
    vsq = vx**2 + vy**2
    W = 1 / np.sqrt(1 - vsq)
    KE = rho * W * (W-1)

    return KE

def GetKETF(anim, frame):
    """
    Retrieves and computes the kinetic energy density for each frame in a two fluid animation.

    Parameters
    ----------
    anim : object
        animation class containing all user def variables
    frame : Array
        Frame from the animation class containing all user def variables at the time we want

    """
    vx1 = frame[anim.variables.index("vx1\n"), 4:-4, 4:-4, 0]
    vy1 = frame[anim.variables.index("vy1\n"), 4:-4, 4:-4, 0]
    vx2 = frame[anim.variables.index("vx2\n"), 4:-4, 4:-4, 0]
    vy2 = frame[anim.variables.index("vy2\n"), 4:-4, 4:-4, 0]
    rho1 = frame[anim.variables.index("rho1\n"), 4:-4, 4:-4, 0]
    rho2 = frame[anim.variables.index("rho2\n"), 4:-4, 4:-4, 0]
    rho = rho1 + rho2
    vx = (vx1 + vx2) / 2
    vy = (vy1 + vy2) / 2
    vsq = vx**2 + vy**2
    W = 1 / np.sqrt(1 - vsq)
    KE = rho * W * (W-1)

    return KE


if __name__ == '__main__':

    #  Model Comparison

    if not 'Ideal' in locals():
        Ideal = Anim('Ideal/HighRes/Data/TimeSeries/UserDef/')
        idealT = Ideal.t.index(min(Ideal.t, key=lambda x : abs(x-3.0)))
        Nideal = Ideal.final.c['nx'] // 2
        KESpecIdeal = getPowerSpectrumSq(Ideal.final, GetKESF(Ideal, Ideal.frame[idealT]))

    if not 'Resistive' in locals():
        Resistive = Anim('Resistive/Sigma10/HighRes/Data/TimeSeries/UserDef/')
        resistiveT = Resistive.t.index(min(Resistive.t, key=lambda x : abs(x-3.0)))
        Nresistive = Resistive.final.c['nx'] // 2
        KESpecResistive = getPowerSpectrumSq(Resistive.final, GetKESF(Resistive, Resistive.frame[resistiveT]))

    if not 'TwoFluid' in locals():
        TwoFluid = Anim('TwoFluid/Sigma10Mu100/HighRes/Data/TimeSeries/UserDef/')
        twoFluidT = TwoFluid.t.index(min(TwoFluid.t, key=lambda x : abs(x-3.0)))
        NtwoFluid = TwoFluid.final.c['nx'] // 2
        KESpecTwoFluid = getPowerSpectrumSq(TwoFluid.final, GetKETF(TwoFluid, TwoFluid.frame[twoFluidT]))


  ### Model Power Spectrum

    fig, axs = plt.subplots(1, 1, sharex=True)
    fig.set_size_inches(6,3)
    fig.tight_layout()

    # Kinetic energy density power
    axs.loglog(np.arange(1, Nideal+1), np.arange(1, Nideal+1)*KESpecIdeal, label=r'$Single \ Fluid \ Ideal$')
    axs.loglog(np.arange(1, Nresistive+1), np.arange(1, Nresistive+1)*KESpecResistive, label=r'$Single \ Fluid \ Resistive$')
    axs.loglog(np.arange(1, NtwoFluid+1), np.arange(1, NtwoFluid+1)*KESpecTwoFluid, label=r'$Two \ Fluid \ Resistive$')
    axs.set_ylabel(r"$k|P_{T}(k)|^2$", {'fontsize':'large'})
    axs.set_xlabel(r'$k$')
    axs.loglog([3, 94.868], [7*10**-2, 7*10**(-2 - 1.5*5/3)], 'k--')
    axs.annotate(r'$k^{-5/3}$', xy=(40, 0.01), fontsize=15)
    axs.set_xlim([1, Nideal])
    axs.legend(loc='lower left')


    plt.savefig('Figures/KineticEnergyPowerSpectrum.eps', format='eps', dpi=1200, bbox_inches='tight')
    plt.show()
