from MMTK import *
from MMTK.Environment import PathIntegrals
from MMTK.Trajectory import Trajectory
from sys import argv
from Scientific import N
from Scientific.Statistics import mean, standardDeviation
from math import sqrt
from numpy import zeros
#from mbpol import mbpolForceField

universe = InfiniteUniverse()
# nanometers

traj = argv[1]

#temperature = 10.0
temperature = float(traj[traj.find("T")+1:traj.find("P")])
P = int(traj[traj.find("P")+1:traj.find("R")])
#P = 4
nmolecules = int(traj[traj.find("N")+1:traj.find("H")])
lattice_spacing = float(traj[traj.find("R")+1:traj.find("F")])
string = str(traj[traj.find("E")+1:traj.find("V")])

if (string == "P"):
    mbpol_test=True
else:
    mbpol_test=False

universe.addObject(PathIntegrals(temperature))

center = zeros( (nmolecules,3) , float)

for i in range(nmolecules):
	universe.addObject(Molecule('water', position = Vector(i*lattice_spacing, 0., 0.)))
        center[i][0] = i*lattice_spacing

for atom in universe.atomList():
	atom.setNumberOfBeads(P)

natoms = len(universe.atomList())

#universe.setForceField(mbpolForceField(universe))

#This is the conversion factor to Units of K
Kper1overcm=11604.505/8065.54445
conv=Kper1overcm/1.196e-2
#print 1./(Units.k_B*0.37*Units.K)

#print traj
#rotskipval=float(argv[2])
#mbpol_test=True
if mbpol_test:
    trajectory = Trajectory(universe, traj)
else:
    trajectory = Trajectory(None, traj)
print 'test'
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
