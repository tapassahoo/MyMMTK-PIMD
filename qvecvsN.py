import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from numpy import zeros,sqrt,mean
import matplotlib.patches as mpatches

valuesx = argv[1]
valuesy = argv[2]
valuesz = argv[3]
numberofsteps = int(valuesx[valuesx.find("N")+1:valuesx.find("S")])

timestepx = zeros( numberofsteps , float )
timestepy = zeros( numberofsteps , float )
timestepz = zeros( numberofsteps , float )

valuesx = open(valuesx, "r")
valuesy = open(valuesy, "r")
valuesz = open(valuesz, "r")

for i in range(numberofsteps):

    valuex = float(valuesx.readline())
    valuey = float(valuesy.readline())
    valuez = float(valuesz.readline())

    timestepx[i] = valuex
    timestepy[i] = valuey
    timestepz[i] = valuez

valuesx.close()
valuesy.close()
valuesz.close()

time = zeros( numberofsteps , float  )

for i in range(numberofsteps):

    time[i] = float(i)/float(numberofsteps)

plt.plot(time, timestepx, 'r-')
plt.plot(time, timestepy, 'b-')
plt.plot(time, timestepz, 'y-')
plt.xlabel('time (scaled)')
plt.ylabel('x,y,z')
plt.title('r vs. time')
plt.show()

