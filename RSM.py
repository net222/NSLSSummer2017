import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import griddata
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from matplotlib.colors import LogNorm

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

def load_RSM_scan(path):
    h = []
    l = []
    I = []
    f = open(path,'r')
    line = f.readline()
    while(line):
        words = line.split(' ')
        if words:
            h.append(float(words[0]))
            l.append(float(words[1]))
            I.append(float(words[2])+1)
        line = f.readline()
    return h, l, I

class RSM_scan(object):
    def __init__(self, path):
        self.import_data(path)
    def import_data(self, path):
        self.h, self.l, self.I = load_RSM_scan(path)

path = r'/Users/Rui/Dropbox/Research/STONY BROOK/Matt Dawber/BTO_growthrate_Project/Raw Data/XRD/J1123161' \
       r'/J1123161_RSM_002_vs_s.dat'
#path = r'/Users/Rui/Dropbox/Research/STONY BROOK/Matt Dawber/PTO_BTO Project/RAW DATA/Anya.dat'
#J0420161 = RSM_scan(path)

h, l, I = load_RSM_scan(path)
z = np.array(I)
x = sorted(list(set(h)))
y = sorted(list(set(l)))
z = z.reshape(len(x),len(y)).T

#levls =np.logspace(np.log(z.min()),np.log(z.max()),30)
width = 10
levls = np.linspace(1,10,width)
levls = np.concatenate((levls[:-1],np.linspace(10,100,width)))
levls = np.concatenate((levls[:-1],np.linspace(100,1000,width)))
levls = np.concatenate((levls[:-1],np.linspace(1000,10000,width)))
levls = np.concatenate((levls[:-1],np.linspace(10000,100000,width)))
levls = np.concatenate((levls[:-1],np.linspace(100000,1000000,width)))
#levls = np.concatenate((levls[:-1],np.linspace(1000000,10000000,width)))
print(z.max())
print(z.shape)


plt.contourf(z,extent = [min(x),max(x),min(y),max(y)],aspect=(max(x)-min(x))/(max(y)-min(y))
           ,cmap=plt.cm.jet,norm = LogNorm(),levels=levls)
#plt.xlim(-0.08,0.08)
#plt.ylim(1.75,2.02)
cbar=plt.colorbar()
cbar.set_ticks([10**0,10**1,10**2,10**3,10**4,10**5,10**6,10**7])

plt.contour(z, extent=[min(x), max(x), min(y), max(y)], aspect=(max(x) - min(x)) / (max(y) - min(y))
             , levels=levls,colors='k',linestyle='.',linewidths=0.2)


#plt.imshow(z, extent=[min(x), max(x), min(y), max(y)], aspect=(max(x) - min(x)) / (max(y) - min(y))
             #, norm=LogNorm(vmin=z.min(), vmax=z.max()), interpolation='nearest')

plt.ylabel("l,l [0 0 1]")
plt.xlabel("h,h [1 0 0]")
plt.title('J1123161_Pure BTO grown at 450C_vs_002 scan')
#plt.tight_layout()
plt.savefig('/Users/Rui/Dropbox/Research/STONY BROOK/Matt Dawber/BTO_growthrate_Project/Raw Data/XRD/J1123161/J1123161_RSM_002_vs_s.pdf')

#plt.autoscale(False)
#plt.plot(xy[:, 1], xy[:, 0], 'ro')


plt.show()










