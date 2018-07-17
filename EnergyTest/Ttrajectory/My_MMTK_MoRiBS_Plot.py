import matplotlib.patches as mpatches 
import matplotlib.pyplot as plt 
from sys import argv                
import numpy as np
from numpy import zeros,sqrt,mean 

erotvalues = np.loadtxt(argv[1])
plt.plot(erotvalues[:,0], erotvalues[:,1], 'r-')
plt.errorbar(erotvalues[:,2], erotvalues[:,3], yerr = erotvalues[:,4])
plt.xlabel('Temperature (K)')
plt.ylabel('Rotational Energy (K)')
plt.axis([0,60,0,80])
plt.title('My MMTK vs. MoRiBS (Rot Only)')
plt.text(5, 60, r'None --> MoRiBS 150 BLOCKS 2000 PASESS')
plt.text(5, 68, r'Blue --> MoRiBS 15 BLOCKS 200 PASSES')
plt.text(5, 76, r'Red --> MMTK 2 000 000 Steps')
plt.show()

