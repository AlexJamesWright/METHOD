"""
Script used for the convergence plot (Figure 5) in Wright & Hawke (2018)

The numbers in the interactive plot argument refer to the resolution of the
domain use for the Field Loop Advection test. 
"""


import numpy as np
from matplotlib import pyplot as plt
from interactivePlot import InteractivePlot
from matplotlib import cm
from matplotlib import ticker


if 'n1212i' not in locals():
    n1212i = InteractivePlot('./1212/Initial/')
if 'n1212f' not in locals():
    n1212f = InteractivePlot('./1212/Final/')
if 'n2424i' not in locals():
    n2424i = InteractivePlot('./2424/Initial/')
if 'n2424f' not in locals():
    n2424f = InteractivePlot('./2424/Final/')
if 'n4848i' not in locals():
    n4848i = InteractivePlot('./4848/Initial/')
if 'n4848f' not in locals():
    n4848f = InteractivePlot('./4848/Final/')
if 'n9696i' not in locals():
    n9696i = InteractivePlot('./9696/Initial/')
if 'n9696f' not in locals():
    n9696f = InteractivePlot('./9696/Final/')
if 'n192192i' not in locals():
    n192192i = InteractivePlot('./192192/Initial/')
if 'n192192f' not in locals():
    n192192f = InteractivePlot('./192192/Final/')


# Convergence
errors=[]
# Only compare the middle region of the domain where fields are non-zero
errors.append(np.sum(np.abs(n1212i.prims[7, 3:-3, 3:-3, 0]-n1212f.prims[7, 3:-3, 3:-3, 0]))/6**2)
errors.append(np.sum(np.abs(n2424i.prims[7, 7:-7, 7:-7, 0]-n2424f.prims[7, 7:-7, 7:-7, 0]))/14**2)
errors.append(np.sum(np.abs(n4848i.prims[7, 15:-15, 15:-15, 0]-n4848f.prims[7, 15:-15, 15:-15, 0]))/30**2)
errors.append(np.sum(np.abs(n9696i.prims[7, 30:-30, 30:-30, 0]-n9696f.prims[7, 30:-30, 30:-30, 0]))/60**2)
errors.append(np.sum(np.abs(n192192i.prims[7, 60:-60, 60:-60, 0]-n192192f.prims[7, 60:-60, 60:-60, 0]))/120**2)


errors = np.array(errors)
dx = np.array([n1212i.c['dx'], n2424i.c['dx'], n4848i.c['dx'], n9696i.c['dx'], n192192i.c['dx']])

fit = np.polyfit(np.log(dx), np.log(errors), 1)
fitXs = np.array([4*10**-3, 0.1])
fitYs = np.power(fitXs, fit[0]) * np.exp(fit[1])


plt.figure()
plt.loglog(dx, errors, '.')
plt.loglog(fitXs, fitYs, 'g--', label=r"$Gradient \approx {:.2f}$".format(fit[0]))
plt.ylim([2*10**-7, 4*10**-4])
plt.xlim([4*10**-3, 0.1])
plt.ylabel(r'$L_1 |B_z^{intial} - B_z^{final}|$', {'fontsize':'large'})
plt.xlabel(r'$\Delta x$')
plt.legend()
plt.savefig('FieldLoopAdvectionConvergence.eps', format='eps', dpi=1200)
plt.show()




# Initial and final plots

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b=int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

fig, axs = plt.subplots(1, 2)
fig.set_figheight(10)
fig.set_figwidth(15)
axs[0].imshow(n128128i.prims[7, 4:-4, 4:-4, 0], extent=(-0.5, 0.5, -0.5, 0.5), cmap=cm.viridis_r, vmin=0.0000, vmax=0.0003, aspect='equal')
cax = axs[1].imshow(n128128f.prims[7, 4:-4, 4:-4, 0], extent=(-0.5, 0.5, -0.5, 0.5), cmap=cm.viridis_r, vmin=0.0000, vmax=0.0003, aspect='equal')
plt.subplots_adjust(wspace=0.05, right=0.8)
cbar_ax = fig.add_axes([0.73, 0.318, 0.015, 0.39])
cbar = fig.colorbar(cax, cax=cbar_ax, format=ticker.FuncFormatter(fmt))
cbar.ax.set_ylabel(r'$B_z$')
axs[0].set_xlabel(r'$x$')
axs[0].set_ylabel(r'$y$')
axs[1].set_xlabel(r'$x$')
axs[1].set_ylabel(r'$y$')
axs[0].set_ylim([-0.5, 0.5])
axs[1].set_ylim([-0.5, 0.5])
axs[0].set_xlim([-0.5, 0.5])
axs[1].set_xlim([-0.5, 0.5])
plt.subplots_adjust(wspace=0.2, right=0.7)
plt.savefig('FieldLoopInitialFinal.eps', format='eps', dpi=270)
plt.show()
