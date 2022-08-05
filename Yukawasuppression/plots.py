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
plt.plot(spectra[:,0],spectra[:,1],label='GR')
plt.plot(spectra[:,0],spectra[:,2],label='nDGP Screening')
plt.xscale('log')
plt.xlabel('wavenumber $k$ $[h\mathrm{Mpc}^{-1}]$')
plt.ylabel('$P_{nl}$')
plt.xlim(0.01,10.0)
plt.yscale('log')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

