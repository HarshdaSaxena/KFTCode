import matplotlib.pyplot as plt
import numpy as np


plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title


spectra = np.loadtxt('temp.d')

plt.figure()
#plt.plot(spectra[:,0], spectra[:,2]/spectra[:,1] - 1,label='nDGP')
plt.plot(spectra[:,0],spectra[:,3], label = '$G_{eff}$ of nDGP')
plt.xscale('log')
plt.xlabel('wavenumber $k$ $[h\mathrm{Mpc}^{-1}]$')
plt.ylabel('$G_{eff}/G_N$')
#plt.ylabel('$P_{nl}/P_{GR} - 1$')
plt.xlim(0.1,10.0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

