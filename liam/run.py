from MMTK import *
from MMTK.ForceFields import Amber99ForceField
from mbpol import mbpolForceField
from MMTK.ForceFields.ForceField import CompoundForceField
from MMTK.Trajectory import Trajectory, TrajectoryOutput, LogOutput, StandardLogOutput
import numpy as np

universe = InfiniteUniverse()
sep_dist = 4.0*Units.Ang
pos1 = Vector(0.0,0.0,0.0)
pos2 = Vector(0.0,sep_dist,0.0)
universe.addObject(Molecule('water', position=pos1))
universe.addObject(Molecule('water', position=pos2))

universe.setForceField(mbpolForceField(universe))

print universe.energy()


