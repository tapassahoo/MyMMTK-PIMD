from MMTK import *

from MMTK import Features
from MMTK.ForceFields.ForceField import CompoundForceField
#from MMTK.ForceFields.ForceFieldTest import gradientTest
#from MMTK.Minimization import ConjugateGradientMinimizer
from MMTK.Environment import PathIntegrals
#from MMTK.NormalModes import VibrationalModes
from MMTK.Trajectory import Trajectory, TrajectoryOutput, \
                            RestartTrajectoryOutput, StandardLogOutput, \
                            trajectoryInfo
from sys import argv,exit
from numpy import zeros,cos,sin,sqrt,pi, dot, asarray, sign, arctan

#CHANGE INTEGRATOR HERE AND FURTHER BELOW WHERE THE INTEGRATOR IS CALLED

from H20Rigid3DRotor_PINormalModeIntegrator import Rigid3DRotor_PINormalModeIntegrator, Rigid3DRotor_PILangevinNormalModeIntegrator

from RotOnlyPINormalModeIntegrator import RotOnly3D_PINormalModeIntegrator, RotOnly3D_PILangevinNormalModeIntegrator

#CHANGE POTENTIAL (FORCEFIELD) HERE AND FURTHER BELOW WHERE THE POTENTIAL IS CALLED

from LJFF import LJForceField
from NoPotFF import NoPotForceField
from dipoleFF import dipoleForceField
from HarmonicWellFF import HarmonicWellForceField
from electricfieldFF import electricForceField

#####################################
## set the number of quantum beads ##
#####################################

label = "Dip"
outdir = "/work/l2mcgrat/trajectoryfiles/"

# Parameters
rhoname = argv[1]
erotname = argv[2]
esqname = argv[3]

Rot_Step = argv[4] #Rot Step
Rot_Skip = argv[5] #Rot Skip Step

temperature = float(rhoname[rhoname.find("T")+1:rhoname.find("t")])*Units.K
P = int(rhoname[rhoname.find("t")+1:rhoname.find(".rho")])
numsteps = 1000

lattice_spacing = 6.0*Units.Ang

ndens = 23588101
denrho = zeros(ndens,float)
denerot = zeros(ndens,float)
denesq = zeros(ndens,float)

universe = InfiniteUniverse()
# nanometers

universe.addObject(PathIntegrals(temperature))

nmolecules = 2

center = zeros( (nmolecules,3) , float)

## They will initially be aligned along the x-axis

for i in range(nmolecules):
	universe.addObject(Molecule('water', position = Vector(i*lattice_spacing, 0., 0.)))
        center[i][0] = i*lattice_spacing

for atom in universe.atomList():
	atom.setNumberOfBeads(P)

# print "ATOMS"
# print  universe.atomList()

natoms = len(universe.atomList())

ff=[]

##################################################
############## COULOMB/LJ POTENTIAL ##############
##################################################

#for im in range(nmolecules):
#        for ia in range(len(universe.objectList()[im].atomList())):
#                for jm in range(im+1,nmolecules):
#                        for ja in range(len(universe.objectList()[jm].atomList())):
# 				ff.append(CoulombForceField(universe.objectList()[im].atomList()[ia],universe.objectList()[jm].atomList()[ja],ia,ja))
#                                ff.append(LJForceField(universe.objectList()[im].atomList()[ia],universe.objectList()[jm].atomList()[ja],ia,ja))

##################################################
############## DIPOLE POTENTIAL ##################
##################################################

#for i in range(nmolecules):
#       for j in range(i+1,nmolecules):
#               ff.append(dipoleForceField(universe.objectList()[i].atomList()[0],universe.objectList()[i].atomList()[1],universe.objectList()[i].atomList()[2],
#                   universe.objectList()[j].atomList()[0],universe.objectList()[j].atomList()[1],universe.objectList()[j].atomList()[2]))

##################################################
############### ELECTRIC FIELD ###################
##################################################

for i in range(nmolecules):
        ff.append(electricForceField(universe.objectList()[i].atomList()[0], universe.objectList()[i].atomList()[1], universe.objectList()[i].atomList()[2]))

##################################################
############## HARMONIC WELL #####################
##################################################
#
# for i in range(nmolecules):
#        ff.append(HarmonicWellForceField(universe.objectList()[i].atomList()[2], center[i][0], center[i][1], center[i][2], float(31.273e4)))

##################################################
############### EMPTY POTENTIAL ##################
##################################################

# for i in range(nmolecules):
#        for j in range(i+1,nmolecules):
#                ff.append(NoPotForceField(universe.objectList()[i].atomList()[ia],universe.objectList()[j].atomList()[ja],ia,ja))

universe.setForceField(CompoundForceField(*ff))

# gradientTest(universe)
# print universe.energy()
# raise()
# print "AFTER SIMULATION"
# universe.setForceField(HarmonicDistanceRestraint(universe[0].H1,universe[0].F1,lattice_spacing, k))

universe.writeToFile("u.pdb")
#universe.environmentObjectList(Environment.PathIntegrals)[0].include_spring_terms = False
#universe._changed(True)
universe.initializeVelocitiesToTemperature(temperature)

rhofile=open(rhoname,"r")
erotfile=open(erotname, "r")
esqfile=open(esqname, "r")

for i in range(ndens):
        denrho[i]=float(rhofile.readline())
        denerot[i]=float(erotfile.readline())*1.1963e-2 #cm-1 to kJ/mol [MMTK Units of Energy]
        denesq[i]=float(esqfile.readline())*(1.1963e-2)*(1.1963e-2)

rhofile.close()
erotfile.close()
esqfile.close()

dt = 1.0*Units.fs
#print "Using a timestep of ", dt*1000., " fs"

# Initialize velocities
universe.initializeVelocitiesToTemperature(temperature)

# USE THE FRICTION PARAMETER FROM BEFORE
friction = 0.0

#integrator = Rigid3DRotor_PILangevinNormalModeIntegrator(universe, delta_t = dt, centroid_friction = friction, denrho = denrho, denerot = denerot, denesq = denesq, rotstep = float(Rot_Step), rotskipstep=float(Rot_Skip))

integrator = RotOnly3D_PILangevinNormalModeIntegrator(universe, delta_t = dt, centroid_friction = friction, denrho = denrho, denerot = denerot, denesq = denesq, rotstep = float(Rot_Step), rotskipstep=float(Rot_Skip))

integrator(steps = 30, actions = [ TrajectoryOutput(None,('configuration','time'), 0, None, 100)] )  # relates to the default_options = {'first_step': 0...} section of the main code.

RunSteps = int(numsteps)*Units.fs/dt
SkipSteps = 1.0*Units.fs/dt

trajectory = Trajectory(universe, outdir+"N"+str(nmolecules)+"H20T"+str(temperature)+"P"+str(P)+"R"+str(lattice_spacing)+"FilEDVersion"+str(numsteps)+"Steps"+label+".nc", "w", "A simple test case")

Nblocks = 1

############################## BEGIN ROTATION SIMULATION ##############################


# RUN PIMD WITH PIMC ROTATION INCLUDED
print "We're going to run the Langevin integrator for ", RunSteps/SkipSteps, "independent steps of PIMD"
integrator(steps=RunSteps,
           # Remove global translation every 50 steps.
	   actions = [
		   TrajectoryOutput(trajectory, ("time", "thermodynamic", "energy",
						 "configuration", "auxiliary"),
                                    0, None, SkipSteps)])

raise()
##############################              BEGIN ANALYSIS           ##############################
#gradientTest(universe)

cosfile1 = open("cosfile1-"+str(nmolecules)+"H20T"+str(temperature)+"P"+str(P)+"FileBVersioN"+str(numsteps)+"Steps"+label,"w")
cosfile2 = open("cosfile2-"+str(nmolecules)+"H20T"+str(temperature)+"P"+str(P)+"FileBVersioN"+str(numsteps)+"Steps"+label,"w")

for step in trajectory:
    universe.setConfiguration(step['configuration'])
    for i in range(P):
        for a in range(len(universe.objectList()) - 1):         # loops over all molecules
            xm1a1=universe.objectList()[a].atomList()[0].beadPositions()[i]
            xm1a2=universe.objectList()[a].atomList()[1].beadPositions()[i]
            xm1a3=universe.objectList()[a].atomList()[2].beadPositions()[i]
            xm2a1=universe.objectList()[a+1].atomList()[0].beadPositions()[i]
            xm2a2=universe.objectList()[a+1].atomList()[1].beadPositions()[i]
            xm2a3=universe.objectList()[a+1].atomList()[2].beadPositions()[i]

            mo1=universe.objectList()[a].atomList()[0].mass()
            mo2=universe.objectList()[a].atomList()[1].mass()
            mo3=universe.objectList()[a].atomList()[2].mass()

            cm1 = (xm1a1*mo1 + xm1a2*mo2 + xm1a3*mo3) / (mo1 + mo2 + mo3)
            cm2 = (xm2a1*mo1 + xm2a2*mo2 + xm2a3*mo3) / (mo1 + mo2 + mo3)

            rcom1 = xm1a3 - cm1
            rcom2 = xm2a3 - cm2

            norm1 = sqrt(dot(rcom1,rcom1))
            norm2 = sqrt(dot(rcom2,rcom2))

            rcom1 /= norm1
            rcom2 /= norm2

            cost1 = rcom1[0]
            cost2 = rcom2[0]

            cosfile1.write(str(cost1)+"\n")
            cosfile2.write(str(cost2)+"\n")

###################################################################################################
##############################       AUTOCORRELATION FUNCTIONS       ##############################

T = len(trajectory) #numbsteps
x = zeros( len(trajectory) , float )
counter = -1

#  we are trying to find C(tau) for the x position of the Oxygen atom on the water molecule
#  for the average bead position of the oxygen in each bead.

for step in trajectory:
    universe.setConfiguration(step['configuration'])
    counter += 1
    x[counter] = universe.objectList()[1].atomList()[2].position()[0]

C1file=open("C1file-"+str(nmolecules)+"-P-"+str(P)+"-"+label,"w")
C2file=open("C2file-"+str(nmolecules)+"-P-"+str(P)+"-"+label,"w")

for tau in range(T):
    C1 = 0.0
    C2 = 0.0

    for t in range(T - tau):
        C1 += ( x[t]*x[t+tau] - mean*mean ) / ( ( T + 1 - tau )*( variance ) )
        C2 += ( x[t] - mean )*( x[t+tau] - mean ) / ( ( T + 1 - tau )*( variance ) )

    if tau <= 3*T/4:
        C1file.write(str(C1)+"\n")
        C2file.write(str(C2)+"\n")

###################################################################################################

npoints = len(trajectory)
universe = trajectory.universe
natoms = universe.numberOfAtoms()
time = trajectory.time
np = universe.numberOfPoints()
P = np / natoms
#gradientTest(universe)
trajectory.close()

afile1.close()
afile2.close()
Xfile.close()
C1file.close()
C2file.close()

