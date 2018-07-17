/* C routines for dipoleFF.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
#include <math.h>

/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
dipole_evaluator(PyFFEnergyTermObject *self,
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
                              gradients, second derivatives).
     */
{
  vector3 *coordinates = (vector3 *)input->coordinates->data;
  vector3 *g;
  int atom_index1 = (int)self->param[0];  /* atom index */
  int atom_index2 = (int)self->param[1];  /* atom index */
  int atom_index3 = (int)self->param[2];  /* atom index */
  int atom_index4 = (int)self->param[3];  /* atom index */
  int atom_index5 = (int)self->param[4];  /* atom index */
  int atom_index6 = (int)self->param[5];  /* atom index */

/*
  printf("%i %i %i %i %i %i \n",atom_index1,atom_index2,atom_index3,atom_index4,atom_index5,atom_index6);
*/

  double x1 = coordinates[atom_index1][0];
  double y1 = coordinates[atom_index1][1];
  double z1 = coordinates[atom_index1][2];

  double x2 = coordinates[atom_index2][0];
  double y2 = coordinates[atom_index2][1];
  double z2 = coordinates[atom_index2][2];

  double x3 = coordinates[atom_index3][0];
  double y3 = coordinates[atom_index3][1];
  double z3 = coordinates[atom_index3][2];

  double x4 = coordinates[atom_index4][0];
  double y4 = coordinates[atom_index4][1];
  double z4 = coordinates[atom_index4][2];

  double x5 = coordinates[atom_index5][0];
  double y5 = coordinates[atom_index5][1];
  double z5 = coordinates[atom_index5][2];

  double x6 = coordinates[atom_index6][0];
  double y6 = coordinates[atom_index6][1];
  double z6 = coordinates[atom_index6][2];

  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term. */

  // unit convertion to MMTK units
  double massO = 15.9993047236;
  double massH = 1.00797597651;
  double e, einK, bondlength, r, charge;
  double massOpercent, massHpercent, mu1dotmu2, mu1dotr, mu2dotr;
  double gr[18];
  double mu1[3];
  double mu2[3];
  double cmvec[3];
  double check;

  bondlength = 0.05882523055; //Bond Length (nm)
  charge = (1.8546 / (bondlength*1.0e-9))*3.33564e-30; //Coulomb

  double q = charge;

  mu1[0] = x3 - (x2 + x1)/2;
  mu1[1] = y3 - (y2 + y1)/2;
  mu1[2] = z3 - (z2 + z1)/2;

  mu2[0] = x6 - (x5 + x4)/2;
  mu2[1] = y6 - (y5 + y4)/2;
  mu2[2] = z6 - (z5 + z4)/2;

  massOpercent = massO/(massO + 2*massH);
  massHpercent = massH/(massO + 2*massH);

  cmvec[0] = (massOpercent*(x6 - x3)) + massHpercent*((x5 + x4) - (x2 + x1));
  cmvec[1] = (massOpercent*(y6 - y3)) + massHpercent*((y5 + y4) - (y2 + y1));
  cmvec[2] = (massOpercent*(z6 - z3)) + massHpercent*((z5 + z4) - (z2 + z1));

  r = sqrt(cmvec[0]*cmvec[0] + cmvec[1]*cmvec[1] + cmvec[2]*cmvec[2]);

  mu1dotmu2 = mu1[0]*mu2[0] + mu1[1]*mu2[1] + mu1[2]*mu2[2];
  mu1dotr = mu1[0]*cmvec[0] + mu1[1]*cmvec[1] + mu1[2]*cmvec[2];
  mu2dotr = mu2[0]*cmvec[0] + mu2[1]*cmvec[1] + mu2[2]*cmvec[2];

  e = (q*q/(4.0*M_PI*r*r*r*8.854187817e-12))*(mu1dotmu2 - 3.0*mu1dotr*mu2dotr/(r*r))*(1.0e9*6.022140857e23)/1000.0;
  check = (mu1dotmu2 - 3.0*mu1dotr*mu2dotr/(r*r));
  einK = e*120.3289;

  printf("%f \n", check);

  double predv=(q*q/(4.0*M_PI*r*r*r*8.854187817e-12))*(1.0e9*6.022140857e23)/1000.0;
  // components

  gr[0] = ((3.0*predv)/(r*r))*cmvec[0]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu2[0] - 3.0*mu2dotr*(cmvec[0] - massHpercent*mu1[0])/(r*r) + 3.0*mu1dotr*mu2[0]*massHpercent/(r*r));
  gr[1] = ((3.0*predv)/(r*r))*cmvec[1]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu2[1] - 3.0*mu2dotr*(cmvec[1] - massHpercent*mu1[1])/(r*r) + 3.0*mu1dotr*mu2[1]*massHpercent/(r*r));
  gr[2] = ((3.0*predv)/(r*r))*cmvec[2]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu2[2] - 3.0*mu2dotr*(cmvec[2] - massHpercent*mu1[2])/(r*r) + 3.0*mu1dotr*mu2[2]*massHpercent/(r*r));

  gr[3] = ((3.0*predv)/(r*r))*cmvec[0]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu2[0] - 3.0*mu2dotr*(cmvec[0] - massHpercent*mu1[0])/(r*r) + 3.0*mu1dotr*mu2[0]*massHpercent/(r*r));
  gr[4] = ((3.0*predv)/(r*r))*cmvec[1]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu2[1] - 3.0*mu2dotr*(cmvec[1] - massHpercent*mu1[1])/(r*r) + 3.0*mu1dotr*mu2[1]*massHpercent/(r*r));
  gr[5] = ((3.0*predv)/(r*r))*cmvec[2]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu2[2] - 3.0*mu2dotr*(cmvec[2] - massHpercent*mu1[2])/(r*r) + 3.0*mu1dotr*mu2[2]*massHpercent/(r*r));

  gr[6] = ((3.0*predv)/(r*r))*cmvec[0]*massOpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*((-1.0)*mu2[0] - 3.0*mu2dotr*((-1.0)*cmvec[0] - massOpercent*mu1[0])/(r*r) + 3.0*mu1dotr*mu2[0]*massOpercent/(r*r));
  gr[7] = ((3.0*predv)/(r*r))*cmvec[1]*massOpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*((-1.0)*mu2[1] - 3.0*mu2dotr*((-1.0)*cmvec[1] - massOpercent*mu1[1])/(r*r) + 3.0*mu1dotr*mu2[1]*massOpercent/(r*r));
  gr[8] = ((3.0*predv)/(r*r))*cmvec[2]*massOpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*((-1.0)*mu2[2] - 3.0*mu2dotr*((-1.0)*cmvec[2] - massOpercent*mu1[2])/(r*r) + 3.0*mu1dotr*mu2[2]*massOpercent/(r*r));

  gr[9] = (-1.0)*((3.0*predv)/(r*r))*cmvec[0]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu1[0] - 3.0*mu2dotr*mu1[0]*massHpercent - 3.0*mu1dotr*(cmvec[0] + massHpercent*mu2[0])/(r*r));
  gr[10] = (-1.0)*((3.0*predv)/(r*r))*cmvec[1]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu1[1] - 3.0*mu2dotr*mu1[1]*massHpercent - 3.0*mu1dotr*(cmvec[1] + massHpercent*mu2[1])/(r*r));
  gr[11] = (-1.0)*((3.0*predv)/(r*r))*cmvec[2]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu1[2] - 3.0*mu2dotr*mu1[2]*massHpercent - 3.0*mu1dotr*(cmvec[2] + massHpercent*mu2[2])/(r*r));

  gr[12] = (-1.0)*((3.0*predv)/(r*r))*cmvec[0]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu1[0] - 3.0*mu2dotr*mu1[0]*massHpercent - 3.0*mu1dotr*(cmvec[0] + massHpercent*mu2[0])/(r*r));
  gr[13] = (-1.0)*((3.0*predv)/(r*r))*cmvec[1]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu1[1] - 3.0*mu2dotr*mu1[1]*massHpercent - 3.0*mu1dotr*(cmvec[1] + massHpercent*mu2[1])/(r*r));
  gr[14] = (-1.0)*((3.0*predv)/(r*r))*cmvec[2]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*(mu1[2] - 3.0*mu2dotr*mu1[2]*massHpercent - 3.0*mu1dotr*(cmvec[2] + massHpercent*mu2[2])/(r*r));

  gr[15] = (-1.0)*((3.0*predv)/(r*r))*cmvec[0]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*((-1.0)*mu1[0] - 3.0*mu2dotr*massOpercent*mu1[0]/(r*r) - 3.0*mu1dotr*((-1.0)*cmvec[0] + mu2[0]*massOpercent)/(r*r));
  gr[16] = (-1.0)*((3.0*predv)/(r*r))*cmvec[1]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*((-1.0)*mu1[1] - 3.0*mu2dotr*massOpercent*mu1[1]/(r*r) - 3.0*mu1dotr*((-1.0)*cmvec[1] + mu2[1]*massOpercent)/(r*r));
  gr[17] = (-1.0)*((3.0*predv)/(r*r))*cmvec[2]*massHpercent*(mu1dotmu2 - 5.0*mu1dotr*mu2dotr/(r*r)) + predv*((-1.0)*mu1[2] - 3.0*mu2dotr*massOpercent*mu1[2]/(r*r) - 3.0*mu1dotr*((-1.0)*cmvec[2] + mu2[2]*massOpercent)/(r*r));

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

  g[atom_index4][0]+=gr[9];
  g[atom_index4][1]+=gr[10];
  g[atom_index4][2]+=gr[11];

  g[atom_index5][0]+=gr[12];
  g[atom_index5][1]+=gr[13];
  g[atom_index5][2]+=gr[14];

  g[atom_index6][0]+=gr[15];
  g[atom_index6][1]+=gr[16];
  g[atom_index6][2]+=gr[17];

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
   module, HeHeFF.py. */
static PyObject *
dipoleTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int atom_index1;
  int atom_index2;
  int atom_index3;
  int atom_index4;
  int atom_index5;
  int atom_index6;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!iiiiii",
			&PyUniverseSpec_Type, &self->universe_spec,
			&atom_index1, &atom_index2, &atom_index3, &atom_index4, &atom_index5, &atom_index6))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = dipole_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "dipole";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("dipole");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */
  self->param[0] = (double) atom_index1;
  self->param[1] = (double) atom_index2;
  self->param[2] = (double) atom_index3;
  self->param[3] = (double) atom_index4;
  self->param[4] = (double) atom_index5;
  self->param[5] = (double) atom_index6;
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
  {"dipoleTerm", dipoleTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_dipole(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_dipole", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_dipole");
}
