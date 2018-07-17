# q-TIP4P water potential
# Kevin Bishop

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
# from MMTK_mbpol import mbpolTerm
from MMTK.Units import electrostatic_energy
from mbpol_eval import mbpolEvalEnergy,mbpolEvalEnergyAndGrad
import numpy as np
from MMTK import ParticleVector
from Scientific.Geometry import Vector

class mbpolTerm(EnergyTerm):

    # The __init__ method only remembers parameters. Note that
    # EnergyTerm.__init__ takes care of storing the name and the
    # universe object.
    def __init__(self, universe, atoms):

        self.mbpol_atom_indices = []
        #Loop over the water molecules, place Oxygen index first, followed by H1,H2 for each water
        for i in range(len(universe.objectList())):
            self.mbpol_atom_indices.append(atoms[i*3+2])
            self.mbpol_atom_indices.append(atoms[i*3])
            self.mbpol_atom_indices.append(atoms[i*3+1])

        EnergyTerm.__init__(self, 'mbpol_potential', universe)

    # This method is called for every single energy evaluation, so make
    # it as efficient as possible. The parameters do_gradients and
    # do_force_constants are flags that indicate if gradients and/or
    # force constants are requested.
    def evaluate(self, configuration, do_gradients, do_force_constants):
        mbpol_configuration = []
        for index in self.mbpol_atom_indices:
            for i in range(3):
                mbpol_configuration.append(configuration.array[index][i]*10.)   #convert to angstroms

        results = {}
        if do_gradients:
            E,g = mbpolEvalEnergyAndGrad(np.array(mbpol_configuration),len(mbpol_configuration)/9)
            results['energy'] = E*4.184
            # results['mbpol_potential'] = E
            gradients = ParticleVector(self.universe)
            # print self.mbpol_atom_indices
            for index in self.mbpol_atom_indices:
                loc = self.mbpol_atom_indices.index(index) * 3
                # print loc
                gradients[index] += Vector(g[loc],g[loc+1],g[loc+2])*4.184*10.
            # print g
            # print gradients.array
            results['gradients'] = gradients
        else:
            E,g = mbpolEvalEnergyAndGrad(np.array(mbpol_configuration),len(mbpol_configuration)/9)
            # results['mbpol_potential'] = E
            results['energy'] = E*4.184
            # results['energy'] = 0.
            # gradients = ParticleVector(self.universe)
            # # print self.mbpol_atom_indices
            # for index in self.mbpol_atom_indices:
            #     loc = self.mbpol_atom_indices.index(index) * 3
            #     # print loc
            #     gradients[index] += Vector(g[loc],g[loc+1],g[loc+2])
            # # print g
            # # print gradients.array
            # results['gradients'] = gradients
        return results

        # results = {}
        # d = configuration[self.atom_index]-self.center
        # results['energy'] = 0.5*self.force_constant*(d*d)
        # if do_gradients:
        #     gradients = ParticleVector(self.universe)
        #     gradients[self.atom_index] = self.force_constant*d
        #     results['gradients'] = gradients
        # if do_force_constants:
        #     force_constants = SymmetricPairTensor(self.universe)
        #     force_constants[self.atom_index, self.atom_index] = self.force_constant*delta
        #     results['force_constants'] = force_constants
        # return results


class mbpolForceField(ForceField):

    """
    mbpol potential from Paesani and co-workers
    """


    def __init__(self,universe):
        self.atoms = map(lambda x: self.getAtomParameterIndices([x])[0],universe.atomList())
        self.arguments = (self.atoms)
        ForceField.__init__(self, 'mbpol_potential')

    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")

        f, offsets = self.beadOffsetsAndFactor(self.atoms, global_data)

        return [mbpolTerm(universe,self.atoms + o)
                for o in offsets]

