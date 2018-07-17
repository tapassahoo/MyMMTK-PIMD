from MMTK import *
from MMTK.Environment import PathIntegrals
from MMTK.Trajectory import Trajectory
from sys import argv
from Scientific import N
from Scientific.Statistics import mean, standardDeviation
from math import sqrt

traj = argv[1]
#print traj
#rotskipval=float(argv[2])
trajectory = Trajectory(None, traj)
print 'test'
