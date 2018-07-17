from numpy import zeros,cos,sin,sqrt,pi, dot, asarray, sign, arctan
from sys import argv,exit

nmolecules = 1
temperature = 50.0
P = 4
numsteps = 1000
label = "TDE"

numbers = open("numbers-"+str(nmolecules)+"H20T"+str(temperature)+"P"+str(P)+"FileBVersioN"+str(numsteps)+"Steps"+label,"w")

for i in range(numsteps):
    value1 = -1.0 + 1.0*(2.0*i + 1.0)/numsteps
    numbers.write(str(value1)+"\n")
    value2 = -1.0 + 2.0*(2.0*i + 1.0)/numsteps
    numbers.write(str(value2)+"\n")
    value3 = -1.0 + 4.0*(2.0*i + 1.0)/numsteps
    numbers.write(str(value3)+"\n")
    value4 = -1.0 + 8.0*(2.0*i + 1.0)/numsteps
    numbers.write(str(value4)+"\n")

numbers.close()
