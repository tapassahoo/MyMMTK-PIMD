import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from numpy import zeros,sqrt,mean
import matplotlib.patches as mpatches

values1 = argv[1]
numberofsteps = int(values1[values1.find("N")+1:values1.find("S")])

timestep1 = zeros( numberofsteps , float )
timestep2 = zeros( numberofsteps , float )
timestep3 = zeros( numberofsteps , float )
timestep4 = zeros( numberofsteps , float )

values1 = open(values1, "r")

for i in range(numberofsteps):

    value1 = float(values1.readline())
    value2 = float(values1.readline())
    value3 = float(values1.readline())
    value4 = float(values1.readline())
    timestep1[i] = value1
    timestep2[i] = value2
    timestep3[i] = value3
    timestep4[i] = value4

values1.close()

time = zeros( numberofsteps , float  )

for i in range(numberofsteps):

    time[i] = float(i)/float(numberofsteps)

plt.plot(time, timestep1, 'r-')
plt.plot(time, timestep2, 'b-')
plt.plot(time, timestep3, 'y-')
plt.plot(time, timestep4, 'p-')
plt.xlabel('cos(theta)')
plt.ylabel('frequency')
plt.title('Angle Distribution')
plt.show()

