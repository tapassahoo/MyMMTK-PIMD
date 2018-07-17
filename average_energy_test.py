from sys import argv
from numpy import zeros,sqrt,mean

################# AVERAGE 1 #################
#                                           #
#    Calculates the average of all the      #
#   energy values in the rot.den_eng file   #
#                                           #
#############################################

Z = 0.0
energy = 0.0
energy_squared = 0.0
rhofile = open(argv[1],"r")
energyfile = open(argv[2],"r")


for i in range(23588101):
    density_matrix_value = rhofile.readline()
    Z += float(density_matrix_value)

rhofile.close()
rhofile = open(argv[1],"r")

for i in range(23588101):
    density_matrix_value = rhofile.readline()
    energy_value = energyfile.readline()
    energy += float(energy_value)*float(density_matrix_value)/Z
    energy_squared += float(energy_value)*float(density_matrix_value)*float(energy_value)*float(density_matrix_value)/(Z**2)

rhofile.close()
energyfile.close()

standard_deviation = sqrt(energy_squared-(energy**2))

print energy, "+/-", standard_deviation

