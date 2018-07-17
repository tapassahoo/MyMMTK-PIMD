# Use this script to compile the C module by running
#
#         pydev setup.py build_ext --inplace
#
# and copying the compiled module from the directory 'build'.

from distutils.core import setup, Extension
import os, sys

compile_args = []
prefix_user='/home/l2mcgrat/MMTK/'
include_dirs= [prefix_user+'include/python2.7/',prefix_user+'downloads/path_integrals/Include/',prefix_user+'include']
libraries=['gfortran']
library_dirs=['.']
from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    import numpy.distutils.misc_util
    include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

setup (name = "dipoleForceField",
       version = "1.0",
       description = "dipole forcefield term for MMTK",

       py_modules = ['dipoleFF'],
       ext_modules = [Extension('MMTK_dipole',
                                ['MMTK_dipole.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs,
                                libraries=libraries,
                                library_dirs=library_dirs)]
       )
