import sys 
import os 
import numpy as np
import math 
from scipy import linalg as LA
from scipy.special import factorial as fac
import Wigner

Jmax = int(sys.argv[1])
    
KelvinPerWaveNumber = 0.6950356           # FROM TAPAS SUPPORT.PY
HartreePerKelvin = 315773.2      # Tapas 
DebyePerAu =  0.393456
AuPerWaveNumber = 219474.6305

beta = 0.02                 # 10 Kelvin 
B = 20.561/KelvinPerWaveNumber         
field = float(sys.argv[2])         # in N / C         
dipolemoment = 1.8265*DebyePerAu*HartreePerKelvin # in kelvin

basissize = 0
for j in range(Jmax+1):
        basissize += 2*j+1

H = np.zeros((basissize,basissize), float) 

# The Kinetic Energy

index1 = 0
for j in range(Jmax+1):            # these two loops act as one loop because 
    for m in range(-j,j+1):         # we require both loops to loop through the entire basis
        index2 = 0
        for jp in range(Jmax+1):              # jp stands for j prime (or j') 
            for mp in range(-jp, jp+1):
                if (j == jp) and (m == mp):
                    H[index1,index2] += B*j*(j+1)         # if you want to turn T_rot off then comment this out  
#                    H[index1,index2] += 0
                index2 =  index2+1
        index1 = index1 + 1

# The Potential Energy

index3 = 0

# if we want potential off, don't uncomment prefactor1=2.0*np.sqrt(np.pi/3.0)*field*dipolemoment
prefactor1=0.0

# if we want potential on, comment prefactor1=2.0*np.sqrt(np.pi/3.0)*field*dipolemoment in; 
prefactor1=2.0*np.sqrt(np.pi/3.0)*field*dipolemoment

for j in range (Jmax+1):
    for m in range(-j,j+1):
        index4 = 0
        for jp in range (Jmax+1):
            for mp in range(-jp,jp+1):
                if (m%2==0):
                    prefactor = prefactor1*np.sqrt(((2.0*j+1.0)*3.0*(2.0*jp+1))/(4.0*np.pi))              # determines if m is odd or even because prefactor includes a (-1)^m term.
                else: 
                    prefactor = prefactor1*(-1.0)*np.sqrt(((2.0*j+1.0)*3.0*(2.0*jp+1))/(4.0*np.pi))

                if (Wigner.Wigner3j(j,1,jp,0,0,0)!=0)and (Wigner.Wigner3j(j,1,jp,-m,0,mp)!=0):
                    H[index3,index4] +=prefactor*Wigner.Wigner3j(j,1,jp,0,0,0)*Wigner.Wigner3j(j,1,jp,-m,0,mp)
                index4 = index4 + 1
        index3 = index3 + 1

[eigvalue,eigvector] = LA.eigh(H)        
zdx = eigvalue.argsort()              
eigvalue = eigvalue[zdx]                
eigvector = eigvector[:,zdx]
#print eigvalue/hartree2Kel
#print eigvalue
partition_func = 0.0
upper=0.0
for i in range(basissize):
    partition_func += np.exp(-beta*eigvalue[i])
for i in range(basissize):
    upper += eigvalue[i]*np.exp(-beta*eigvalue[i])

print upper/partition_func 
