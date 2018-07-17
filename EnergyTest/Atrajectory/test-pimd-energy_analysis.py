from MMTK import *
from MMTK.Environment import PathIntegrals
from MMTK.Trajectory import Trajectory
from sys import argv
from Scientific import N
from Scientific.Statistics import mean, standardDeviation
from math import sqrt

#This is the conversion factor to Units of K
Kper1overcm=11604.505/8065.54445
conv=Kper1overcm/1.196e-2
#print 1./(Units.k_B*0.37*Units.K)

traj = argv[1]
#print traj
#rotskipval=float(argv[2])
trajectory = Trajectory(None, traj)
npoints = len(trajectory)
universe = trajectory.universe
natoms = universe.numberOfAtoms()
time = trajectory.time
np = universe.numberOfPoints()
P = np/natoms

#if (rotskipval < 100.0):
#    rotskipratio=1.0
#else:
#    rotskipratio=100.0/rotskipval

# Print averages of the quantu energy estimator
print "Number of Atoms:", natoms
print "Number of Beads:", P
print "Number of Steps:", npoints
#data = trajectory.temperature
#print "Temperature:", mean(data), "+/-",  standardDeviation(data)/sqrt(npoints), " kj/mol"
data = trajectory.quantum_energy_primitive
print "Primitive estimator:", mean(data/Units.k_B), "+/-",  standardDeviation(data/Units.k_B)/sqrt(npoints), " K"
data2 = trajectory.potential_energy
print "Potential estimator:", mean(data2/Units.k_B), "+/-",  standardDeviation(data2/Units.k_B)/sqrt(npoints), " K"
print "Kinetic estimator:", mean(data-data2)/Units.k_B
#data = trajectory.kinetic_energy
#print "Kinetic estimator:", mean(data/Units.k_B), "+/-",  standardDeviation(data/Units.k_B)/sqrt(npoints), " K"
#data = trajectory.spring_energy
#print "Spring estimator:", mean(data/Units.k_B), "+/-",  standardDeviation(data/Units.k_B)/sqrt(npoints), " K"

#data = trajectory.quantum_energy_virial
#print "Virial estimator:", mean(data), "+/-",  standardDeviation(data)/sqrt(npoints), " kj/mol"
#data = trajectory.quantum_energy_centroid_virial
#print "Centroid virial estimator:", mean(data), "+/-",  standardDeviation(data)/sqrt(npoints), " kj/mol"
data3 = trajectory.quantum_energy_rotation
output_rot=open('KErot','w')
for e in data3:
    output_rot.write(str(e)+'\n')
output_rot.close()
print "Rotational Energy estimator:", mean(data3/Units.k_B), "+/-",  standardDeviation(data3/Units.k_B)/sqrt(npoints), " K"
trajectory.close()
