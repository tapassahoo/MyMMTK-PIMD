# Quartic potential with respect to a fixed point in space
from MMTK.ForceFields.ForceField import ForceField
from MMTK_h2_h2o import H2H2OTerm
from numpy import zeros

class H2H2OForceField(ForceField):

    """H2H2O potential with respect to a fixed point in space

    Constructor: H2H2OForceField(|atom|, |center| )

    Arguments:

    |atom| -- an atom object or an integer atom index, specifying the
              atom on which the force field acts

    |center| -- x,y,z coordinates of oxygen, hyd1, hyd2
                in water molecule respectively (9 parameters total)

    """
    def __init__(self, atom1):
        #self.id1 = self.getAtomParameterIds((atom1,))[0]
        self.atom_index = self.getAtomParameterIndices([atom1])[0]
        # From HarmonicOscillatorFF.py

        #self.arguments = (self.id1, beads)     Original Comment out
        self.arguments = (self.atom_index,)
        #self.arguments=(self.id1,)             Previously used

        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'h2h2o')
        # Store the parameters for later use.
 #       self.center = center
       # self.beads=beads

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
	return True
    
    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        #atom1 = self.getAtomParameters(universe, [self.id1], global_data)[0]
        if subset1 is not None:
            s1 = subset1.atomList()
            s2 = subset2.atomList()
            if not ((atom1 in s1 and atom1 in s2)):
                return []


        # INITIALIZE ARRAY for 501 rad points, 181 theta points, 91 chi points
        Nwater=24
        pot=zeros(8251971)
        dvdr=zeros(8251971)
        dvdt=zeros(8251971)
        dvdc=zeros(8251971)
        r_com=zeros(3*Nwater)
        rotmat=zeros(9*Nwater)

        # FILL ARRAYS
        file=open('../pot','r')
        for i in range (8251971):
            dummy=file.readline()
            pot[i]=float(dummy)
        file.close()

        file=open('../gr','r')
        for i in range (8251971):
            dummy=file.readline()
            dvdr[i]=float(dummy)
        file.close()

        file=open('../gt','r')
        for i in range (8251971):
            dummy=file.readline()
            dvdt[i]=float(dummy)
        file.close()

        file=open('../gc','r')
        for i in range (8251971):
            dummy=file.readline()
            dvdc[i]=float(dummy)
        file.close()

        file=open('watercom.24','r')
        for i in range (3*Nwater):
            dummy=file.readline()
            r_com[i]=float(dummy)
        file.close()

        file=open('rotmat.24','r')
        for i in range (9*Nwater):
            dummy=file.readline()
            rotmat[i]=float(dummy)
        file.close()


        f, offsets = self.beadOffsetsAndFactor([self.atom_index], global_data)    
        return f*[H2H2OTerm(universe._spec,
                            pot,dvdr,dvdt,dvdc,r_com, rotmat,
                            self.atom_index + o) for o, in offsets]    
    
    #        return [H2H2OTerm(universe._spec, atom1.index)]
    #        return [H2H2OTerm(universe._spec, atom1.index, self.beads)]

    # This method returns the string that is inserted into the universe
    # descriptions in trajectories. It is the class name followed by
    # the arguments, just what it takes to re-create an equivalent object.
    def description(self):
        return "H2H2OFF.H2H2OForceField" + `self.arguments`
