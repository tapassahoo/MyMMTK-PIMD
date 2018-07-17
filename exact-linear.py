import sys 
import os 
import numpy as np
import math 
from scipy import linalg as LA
from scipy.special import factorial as fac
import Wigner

KelvinPerWaveNumber = 0.6950356
B = 20.561/KelvinPerWaveNumber         # the rotational constant has units of Kelvin now
T = 50.0

lmax = 9

partition_func = 0.0
upper=0.0

for l in range(lmax):
    partition_func += (2l+1)*np.exp(-B*l*(l+1)/T)

for l in range(lmax):
    upper += B*l*(l+1)*(2l+1)*np.exp(-B*l*(l+1)/T)

print upper/partition_func 
