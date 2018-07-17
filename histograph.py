import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from numpy import zeros,sqrt,mean
import matplotlib.patches as mpatches

values1 = argv[1]
values2 = argv[2]
values3 = argv[3]
numberofbins = int(argv[4])
numberofsteps = int(values1[values1.find("N")+1:values1.find("S")])
numberofbeads = int(values1[values1.find("P")+1:values1.find("F")])
filesize = numberofsteps*numberofbeads
bins1 = zeros( numberofbins , float )
bins2 = zeros( numberofbins , float )
bins3 = zeros( numberofbins , float )

values1 = open(values1, "r")

for i in range(filesize):
    value1 = float(values1.readline())
    for j in range(numberofbins):
        if (value1 >= (-1.0 + (2.0*j/numberofbins))) and (value1 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins1[j] += 1.0/filesize

values1.close()

values2 = open(values2, "r")

for i in range(filesize):
    value2 = float(values2.readline())
    for j in range(numberofbins):
        if (value2 >= (-1.0 + (2.0*j/numberofbins))) and (value2 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins2[j] += 1.0/filesize

values2.close()

values3 = open(values3, "r")

for i in range(filesize):
    value3 = float(values3.readline())
    for j in range(numberofbins):
        if (value3 >= (-1.0 + (2.0*j/numberofbins))) and (value3 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins3[j] += 1.0/filesize

values3.close()

costheta1 = zeros( numberofbins , float )
costheta2 = zeros( numberofbins , float )
dihedral = zeros( numberofbins , float )

for i in range(numberofbins):
    costheta1[i] = -1.0 + (2.0*i + 1.0)/numberofbins
    costheta2[i] = -1.0 + (2.0*i + 1.0)/numberofbins
    dihedral[i] = -1.0 + (2.0*i + 1.0)/numberofbins


plt.plot(costheta1, bins1, 'r-')
plt.plot(costheta2, bins2, 'g-')
plt.plot(dihedral, bins3, 'b-')
plt.xlabel('cos(theta)')
plt.ylabel('frequency')
plt.title('Angle Distribution')

red_patch = mpatches.Patch(color='red', label='costheta1')
green_patch = mpatches.Patch(color='blue', label='costheta2')
blue_patch = mpatches.Patch(color='green', label='dihedral angle')

plt.legend(handles=[red_patch, green_patch, blue_patch])
plt.show()

