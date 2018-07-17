import os
from MMTK import *
from mbpol import mbpolForceField
#from qTIP4pFF import HarmonicAngleForceField, HarmonicBondForceField, QuarticBondForceField, ElectrostaticForceField, LennardJonesForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest, virialTest
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
#from MMTK.ForceFields import Amber99ForceField
from MMTK.ForceFields.ForceField import CompoundForceField
from MMTK.Dynamics import TranslationRemover, RotationRemover
from MMTK_PIGSNormalModeIntegrator import PIGSNormalModeIntegrator, PIGSLangevinNormalModeIntegrator
from MMTK_PIGSCartesianIntegrator import PIGSCartesianIntegrator, PIGSLangevinCartesianIntegrator
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.NormalModes import VibrationalModes
from MMTK import Features
from sys import argv
from MMTK.Minimization import SteepestDescentMinimizer
from Scientific import N
from numpy import zeros, array, sqrt, pi, exp, dot, sinh, cosh, power
from Scientific.Statistics import mean, standardDeviation

if len(argv) != 1:
	raise("usgage:"+argv[0]+"\n")

########################### CHOOSE SYSTEM PARAMETERS #########################

#strbeta="03"
#beta_half=0.015/Units.k_B

strbeta="p5"
beta_half=0.25/Units.k_B

dt=0.12*Units.fs                  # FROM STEP 4 : TS
friction=1.0/(0.3*Units.ps)      # FROM STEP 2 : NVE
skip=int(0.055*Units.ps/dt)        # FROM STEP 3 : NVT-SS
nvt_time=int(550.0*Units.ps/dt)
##############################################################################
tau=float(0.0007)
inttau=int(round(float(tau*100000.0)))
print inttau
tau=tau/Units.k_B
P_half=int(beta_half/tau)
P=2*P_half+1
beta=P*tau

print "P= ",P
print "beta ",beta
print "tau ",tau
print "k_B ",Units.k_B

temperature=1./(Units.k_B*beta)

bootstrap_dt=.001*Units.ps
bootstrap_friction = 0.01/bootstrap_dt
universe =InfiniteUniverse()

sep_dist = 4.0*Units.Ang
pos1 = Vector(0.0,0.0,0.0)
pos2 = Vector(0.0,sep_dist,0.0)
label="dim-mbpol-1T-beta"+strbeta+"-tau-"+str(inttau)
universe.addObject(Molecule('spcfw-q', position=pos1))
universe.addObject(Molecule('spcfw-q', position=pos2))

for atom in universe.atomList():
        atom.setNumberOfBeads(P)

universe.addObject(Environment.PathIntegrals(temperature, False))
#BUILD ALL OF THE FORCEFIELD COMPONENTS
#ff = []

#ff.append(mbpolForceField(universe))

#universe.setForceField(CompoundForceField(*ff))
universe.setForceField(mbpolForceField(universe))
#e, g = universe.energyAndGradients()
#print universe.energyTerms()
#print e
#gradientTest(universe)

# HERE WE PERFORM AN ENERGY MINIMIZATION
print "E before steepest descent minimization", universe.energy()
minimizer = SteepestDescentMinimizer(universe, actions=[])
minimizer(convergence = 1.e-8, steps =50)
print "E after steepest descent minimization", universe.energy()

# HERE WE USE NORMAL MODES TO DETERMINE THE TIME-STEP AND FRICTION
# OF OUR SYSTEM
universe.initializeVelocitiesToTemperature(temperature)

integrator = PIGSLangevinNormalModeIntegrator(universe, delta_t=bootstrap_dt,
					      centroid_friction=friction)

integrator(steps=500, actions = [
                               TrajectoryOutput(None,('configuration','time'), 0, None, 100)] )

print "bootstrap_dt",bootstrap_dt

universe.initializeVelocitiesToTemperature(temperature)
print "new dt=",dt, " new friction=", friction

universe.environmentObjectList(Environment.PathIntegrals)[0].include_spring_terms = False
universe._changed(True)

universe.initializeVelocitiesToTemperature(temperature)
#re-create integrator to use calculated friction value
integrator = PIGSLangevinNormalModeIntegrator(universe, delta_t=dt,
                                                centroid_friction = friction)

print "Number of Steps : "+str(nvt_time/skip)


# NOW WE BUILD A TRAJECTORY, AND RUN OUR ACTUAL DYNAMICS!
trajectoryNVE = Trajectory(universe, "/scratch/mdgschmi/"+label+".nc", "w", "PIMD simulation NVE using Langevin Cartesian Integrator")

# NOTE THERE IS A BRIEF INITIALIZATIONS OF 0.05/dt time-steps!

integrator(steps=nvt_time, actions = [
                               TrajectoryOutput(trajectoryNVE,('configuration','energy','time','velocities'), 0, None, skip)] )


trajectoryNVE.close()
os.system("mv /scratch/mdgschmi/"+label+".nc /warehouse/mdgschmi/MBpolDimer/.")
