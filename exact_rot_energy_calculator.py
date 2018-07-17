# This code calculates <E> by performing the sum of E[i]*p[i]/Z up until the 3rd value of J.
# 
# In the asymetric case there appears to be no degeneracy, instead each of the 2J + 1 values
# has their own seperate energy value calculated with the rotational constants (A, B, C).
#
# Information regarding the energy levels can be found from;
# http://www.chemie.unibas.ch/~tulej/Spectroscopy_related_aspects/Lecture31_Spec_Rel_Asp.pdf

import numpy as np

maxJ = 3

# initialization of elements of the sum

p = np.zeros( ((maxJ+1), (2*maxJ+1)) , float )
energylevel = np.zeros( ((maxJ+1), (2*maxJ+1)) , float )

# rotation constants (in cm^-1) which are required for p[i] and E[i]

A = 27.877 
B = 14.512
C = 9.285

h = 6.62606957*(10**-34)
c = 2.998*(10**10)  # cm/s
kB = 1.38064852*(10**-23)

# convert to Joules

A *= h*c
B *= h*c
C *= h*c

# ask for temperature

T = 10.0 
beta = 1 / (kB*T)

# time to calulcate all of the energy levels and their respective probabilities 

p[0][0] = 1.0
p[1][0] = np.exp(-beta*(A + B))
p[1][1] = np.exp(-beta*(A + C))
p[1][2] = np.exp(-beta*(B + C))
p[2][0] = np.exp(-beta*(2*(A + B + C) + 2*np.sqrt((B - C)**2 + (A - C)*(A - B)))) 
p[2][1] = np.exp(-beta*(4*A + B + C))
p[2][2] = np.exp(-beta*(A + 4*B + C))
p[2][3] = np.exp(-beta*(A + B + 4*C))
p[2][4] = np.exp(-beta*(2*(A + B + C) - 2*np.sqrt((B - C)**2 + (A - C)*(A - B))))
p[3][0] = np.exp(-beta*(5*(A + B) + 2*C + 2*np.sqrt(4*((A - B)**2) + (A - C)*(B - C))))
p[3][1] = np.exp(-beta*(5*(A + C) + 2*B + 2*np.sqrt(4*((A - C)**2) - (A - B)*(B - C))))
p[3][2] = np.exp(-beta*(5*(B + C) + 2*A + 2*np.sqrt(4*((B - C)**2) + (A - B)*(A - C))))
p[3][3] = np.exp(-beta*(4*(A + B + C)))
p[3][4] = np.exp(-beta*(5*(A + B) + 2*C - 2*np.sqrt(4*((A - B)**2) + (A - C)*(B - C))))
p[3][5] = np.exp(-beta*(5*(A + C) + 2*B - 2*np.sqrt(4*((A - C)**2) - (A - B)*(B - C))))
p[3][6] = np.exp(-beta*(5*(B + C) + 2*A - 2*np.sqrt(4*((B - C)**2) + (A - B)*(A - C))))

energylevel[0][0] = 0.0
energylevel[1][0] = A + B
energylevel[1][1] = A + C
energylevel[1][2] = B + C
energylevel[2][0] = 2*(A + B + C) + 2*np.sqrt((B - C)**2 + (A - C)*(A - B))
energylevel[2][1] = 4*A + B + C
energylevel[2][2] = A + 4*B + C
energylevel[2][3] = A + B + 4*C
energylevel[2][4] = 2*(A + B + C) - 2*np.sqrt((B - C)**2 + (A - C)*(A - B))
energylevel[3][0] = 5*(A + B) + 2*C + 2*np.sqrt(4*((A - B)**2) + (A - C)*(B - C))
energylevel[3][1] = 5*(A + C) + 2*B + 2*np.sqrt(4*((A - C)**2) - (A - B)*(B - C))
energylevel[3][2] = 5*(B + C) + 2*A + 2*np.sqrt(4*((B - C)**2) + (A - B)*(A - C))
energylevel[3][3] = 4*(A + B + C)
energylevel[3][4] = 5*(A + B) + 2*C - 2*np.sqrt(4*((A - B)**2) + (A - C)*(B - C))
energylevel[3][5] = 5*(A + C) + 2*B - 2*np.sqrt(4*((A - C)**2) - (A - B)*(B - C))
energylevel[3][6] = 5*(B + C) + 2*A - 2*np.sqrt(4*((B - C)**2) + (A - B)*(A - C))

# now we calculate the partition function

Z = 0.0

for i in range(maxJ + 1):
    for j in range(2*maxJ + 1):
        Z += p[i][j]

# now we use this and the energy levels to calculate the rotational energy value

Rotational_Energy = 0.0

for i in range(maxJ + 1):
    for j in range(2*maxJ + 1):
        Rotational_Energy += energylevel[i][j]*p[i][j]/Z

print " Exact Rotational Energy: ", Rotational_Energy, "Joules "
print " Exact Rotational Energy: ", Rotational_Energy*6.022140758e20, "kJ/mol "
print " Exact Rotational Energy: ", Rotational_Energy*6.022140758e23/8.3144621, "K"
