import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title


x,y,function = np.loadtxt('temp.d',delimiter = ',', dtype = float, unpack = True)
func2= function.reshape(200,200)

c = plt.imshow( func2, norm= colors.LogNorm(vmin=func2.min(), vmax= func2.max()), extent=[0,1,0,1])
plt.colorbar(c)


plt.xlabel('scale factor $x$ ')
plt.ylabel('wavenumber ratio $y = alpha/k$')
plt.show()
