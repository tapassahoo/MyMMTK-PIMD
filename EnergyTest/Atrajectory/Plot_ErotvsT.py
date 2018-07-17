import matplotlib.patches as mpatches 
import matplotlib.pyplot as plt 
from sys import argv                
import numpy as np
from numpy import zeros,sqrt,mean 

erotvalues = np.loadtxt(argv[1])
plt.plot(erotvalues[:,0], erotvalues[:,1], 'r-')
plt.errorbar(erotvalues[:,2], erotvalues[:,3], yerr = erotvalues[:,4])
plt.errorbar(erotvalues[:,5], erotvalues[:,6], yerr = erotvalues[:,7])
# plt.errorbar(erotvalues[:,8], erotvalues[:,9], yerr = erotvalues[:,10])
plt.xlabel('Temperature (K)')
plt.ylabel('Rotational Energy (K)')
plt.axis([0,60,0,80])
plt.title('Rotational Energy vs. Temperature')
plt.text(5, 52, r'Red --> MoRiBS Results,')  
plt.text(5, 60, r'Blue --> MMTK 200 000 Steps')
plt.text(5, 68, r'Yellow --> MMTK 2 000 000 Steps')
plt.show()

