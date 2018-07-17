# Quartic potential with respect to a fixed point in space
from MMTK.ForceFields.ForceField import ForceField
from MMTK_electricfield import electricfieldTerm

class electricForceField(ForceField):

    def __init__(self, atom1, atom2, atom3):
							                                                            
        self.atom_index1, self.atom_index2, self.atom_index3 = self.getAtomParameterIndices((atom1,atom2, atom3))
                                                                                                                        
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (self.atom_index1, self.atom_index2, self.atom_index3)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'Electric Field')

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    # The following method is called by the energy evaluation engine
    # to obtain a list of the low-level evaluator objects (the C routines)
    # that handle the calculations.
    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        # The subset evaluation mechanism does not make much sense for
        # this force field, so we just signal an error if someone
        # uses it by accident.
        if subset1 is not None or subset2 is not None:
            raise ValueError, "sorry, no subsets here"
        # Here we pass all the parameters as "simple" data types to
        # the C code that handles energy calculations.
#        return [dipoleTerm(universe._spec,
#                                       self.atom1,self.atom2)]
                                                                                                        
        f, offsets = self.beadOffsetsAndFactor([self.atom_index1, self.atom_index2, self.atom_index3], global_data)
        return [electricfieldTerm(universe._spec,
                            self.atom_index1 + o1,
                            self.atom_index2 + o2, 
                            self.atom_index3 + o3) for o1,o2,o3 in offsets]

    #        return [H2H2OTerm(universe._spec, atom1.index)]
    #        return [H2H2OTerm(universe._spec, atom1.index, self.beads)]

    # This method returns the string that is inserted into the universe
    # descriptions in trajectories. It is the class name followed by
    # the arguments, just what it takes to re-create an equivalent object.
    def description(self):
        return "electricfieldFF.electricForceField" + `self.arguments`
        #return self.__class__.__name__ + `self.arguments`
