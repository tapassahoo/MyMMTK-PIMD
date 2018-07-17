import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from numpy import zeros,sqrt,mean
import matplotlib.patches as mpatches

values1 = argv[1]
numberofbins = int(argv[2])
numberofsteps = int(values1[values1.find("N")+1:values1.find("S")])
numberofbeads = int(values1[values1.find("P")+1:values1.find("F")])
bins1 = zeros( numberofbins , float )
bins2 = zeros( numberofbins , float )
bins3 = zeros( numberofbins , float )
bins4 = zeros( numberofbins , float  )

values1 = open(values1, "r")

for i in range(numberofsteps):
    value1 = float(values1.readline())
    value2 = float(values1.readline())
    value3 = float(values1.readline())
    value4 = float(values1.readline())
    for j in range(numberofbins):
        if (value1 >= (-1.0 + (2.0*j/numberofbins))) and (value1 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins1[j] += 1.0/numberofsteps
        if (value2 >= (-1.0 + (2.0*j/numberofbins))) and (value2 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins2[j] += 1.0/numberofsteps
        if (value3 >= (-1.0 + (2.0*j/numberofbins))) and (value3 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins3[j] += 1.0/numberofsteps
        if (value4 >= (-1.0 + (2.0*j/numberofbins))) and (value4 < ((-1.0 + (2.0*(j+1)/numberofbins)))):
            bins4[j] += 1.0/numberofsteps

values1.close()

costheta1 = zeros( numberofbins , float )
costheta2 = zeros( numberofbins , float )
costheta3 = zeros( numberofbins , float )
costheta4 = zeros( numberofbins , float )

for i in range(numberofbins):
    costheta1[i] = -1.0 + (2.0*i + 1.0)/numberofbins
    costheta2[i] = -1.0 + (2.0*i + 1.0)/numberofbins
    costheta3[i] = -1.0 + (2.0*i + 1.0)/numberofbins
    costheta4[i] = -1.0 + (2.0*i + 1.0)/numberofbins

plt.plot(costheta1, bins1, 'r-')
plt.plot(costheta2, bins2, 'g-')
plt.plot(costheta3, bins3, 'y-')
plt.plot(costheta4, bins4, 'p-')
plt.xlabel('cos(theta)')
plt.ylabel('frequency')
plt.title('Angle Distribution')
plt.show()

