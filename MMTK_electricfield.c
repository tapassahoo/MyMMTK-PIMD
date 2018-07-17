/* C routines for electricfieldFF.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
#include <math.h>

/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */

static void
electricfield_evaluator(PyFFEnergyTermObject *self,
		   PyFFEvaluatorObject *eval,
		   energy_spec *input,
		   energy_data *energy)

/* The four parameters are pointers to structures that are
    defined in MMTK/forcefield.h.

    PyFFEnergyTermObject: All data relevant to this particular
    energy term.

    PyFFEvaluatorObject:  Data referring to the global energy
    evaluation process, e.g. parallelization
    options. Not used here.

    energy_spec:          Input parameters for this routine, i.e.
    atom positions and parallelization parameters.

    energy_data:          Storage for the results (energy terms,
    gradients, second derivatives).  */

{
  vector3 *coordinates = (vector3 *)input->coordinates->data;
  vector3 *g;
  int atom_index1 = (int)self->param[0];  /* H_L atom index */
  int atom_index2 = (int)self->param[1];  /* H_R atom index */
  int atom_index3 = (int)self->param[2];  /* O atom index   */

/*  printf("%i %i %i %i %i %i \n",atom_index1,atom_index2,atom_index3,atom_index4,atom_index5,atom_index6); */

/*  There is N water molecules in this code, each with 3 atoms. This is why there are 3N position vectors.
    I use the positions of these atoms to define the direction of the dipole moment vector.
    The oxygen is used to find the position of the partial negative charge, and the hydrogens are used to find the position of the partial positive charge. */

  double x1 = coordinates[atom_index1][0];
  double y1 = coordinates[atom_index1][1];  /*  First Left Hydrogen  *//*      O1      */
  double z1 = coordinates[atom_index1][2];

  double x2 = coordinates[atom_index2][0];
  double y2 = coordinates[atom_index2][1];  /*  First Right Hydrogen *//*      HL1     */
  double z2 = coordinates[atom_index2][2];

  double x3 = coordinates[atom_index3][0];
  double y3 = coordinates[atom_index3][1];  /*     First Oxygen      *//*      HR1     */
  double z3 = coordinates[atom_index3][2];

  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term. */

  double U, e;
  double costheta;
  double gr[18];

/*  Reference: https://chem.libretexts.org/Core/Physical_and_Theoretical_Chemistry/Physical_Properties_of_Matter/Atomic_and_Molecular_Properties/Intermolecular_Forces/Specific_Interactions/Dipole-Dipole_Interactions */

  double Na = 6.022140857e23; /* units of (1 / mol) we will multiply this by U which is in units of nJ to obtain nJ / mol                 */
  double unitdipole[3];
  double normalizationfactor;
  double dipolemoment = 1.8546*3.33564e-30; // Coulomb*Meters
  double electricfield = 1.0e8;
  double axialfielddirection[3];

/*  We want the first term in front of U to be in units of (kj/mol)*(nm^3) so when it is multiplied by the proceeding                   */
/*  terms (the next is unitless and the one after will have units of nm^-3) it will have units of kj/mol.                              */

  unitdipole[0] = (x1 + x2)/2.0 - x3;  /* Difference between the positively charged base of  */  /*      O                                             */
  unitdipole[1] = (y1 + y2)/2.0 - y3;  /* the ith molecule (the point between the Hydrogens) */  /*     / \     B is the base, mu points from O to B   */
  unitdipole[2] = (z1 + z2)/2.0 - z3;  /* and the Oxygen of the ith molecule.                */  /*    H B H                                           */

  normalizationfactor = sqrt( unitdipole[0]*unitdipole[0] + unitdipole[1]*unitdipole[1] + unitdipole[2]*unitdipole[2] );

  unitdipole[0] = unitdipole[0] / normalizationfactor;
  unitdipole[1] = unitdipole[1] / normalizationfactor;
  unitdipole[2] = unitdipole[2] / normalizationfactor;

  axialfielddirection[0] = 1.0;
  axialfielddirection[1] = 0.0;
  axialfielddirection[2] = 0.0;

  costheta = unitdipole[0]*axialfielddirection[0] + unitdipole[1]*axialfielddirection[1] + unitdipole[2]*axialfielddirection[2];

  U = dipolemoment*electricfield*costheta*Na/1000;
  e = U;

  /* Because its not a coulombic potential, and the only distance dependant term is r which depends on the center of mass position, we end up computing the same result for HL1 and HR1. */
  /* This is why gr[0] = gr[3], gr[1] = gr[4], and gr[2] = gr[5]. Same goes for HL2 and HR2. gr[9] = gr[12], gr[10] = gr[13], gr[11] = gr[14]. */

  gr[0] = 0.0;
  gr[1] = 0.0;
  gr[2] = 0.0;
  gr[3] = 0.0;
  gr[4] = 0.0;
  gr[5] = 0.0;
  gr[6] = 0.0;
  gr[7] = 0.0;
  gr[8] = 0.0;

  energy->energy_terms[self->index] = e;
  /* If only the energy is asked for, stop here. */
  if (energy->gradients == NULL)
    return;

  /* Add the gradient contribution to the global gradient array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased. */

  g = (vector3 *)((PyArrayObject*)energy->gradients)->data;

  g[atom_index1][0]+=gr[0];
  g[atom_index1][1]+=gr[1];
  g[atom_index1][2]+=gr[2];

  g[atom_index2][0]+=gr[3];
  g[atom_index2][1]+=gr[4];
  g[atom_index2][2]+=gr[5];

  g[atom_index3][0]+=gr[6];
  g[atom_index3][1]+=gr[7];
  g[atom_index3][2]+=gr[8];

}

/* A utility function that allocates memory for a copy of a string */
static char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/* The next function is meant to be called from Python. It creates the
   energy term object at the C level and stores all the parameters in
   there in a form that is convient to access for the C routine above.
   This is the routine that is imported into and called by the Python
   module, HeHeFF.py. --> you mean dipoleFF.py ?? */

static PyObject *
electricfieldTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int atom_index1;
  int atom_index2;
  int atom_index3;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;

  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!iii",
			&PyUniverseSpec_Type, &self->universe_spec,
			&atom_index1, &atom_index2, &atom_index3))
    return NULL;

  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */

  Py_INCREF(self->universe_spec);

  /* A pointer to the evaluation routine. */

  self->eval_func = electricfield_evaluator;

  /* The name of the energy term object. */

  self->evaluator_name = "electric field";

  /* The names of the individual energy terms - just one here. */

  self->term_names[0] = allocstring("electric field");

  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;

  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */

  self->param[0] = (double) atom_index1;
  self->param[1] = (double) atom_index2;
  self->param[2] = (double) atom_index3;

//  self->param[1] = (double) atom_index2;
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */

static PyMethodDef functions[] = {
  {"electricfieldTerm", electricfieldTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */

DL_EXPORT(void)
initMMTK_electricfield(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_electricfield", functions);
  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_electricfield");
}

