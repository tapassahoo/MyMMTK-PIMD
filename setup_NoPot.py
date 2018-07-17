# Use this script to compile the C module by running
#
#         pydev setup.py build_ext --inplace
#
# and copying the compiled module from the directory 'build'.

from distutils.core import setup, Extension
import os, sys

compile_args = []
#include_dirs = ['/home/mdgschmi/mmtk_with_pi2_pigs_for_student_new/Include','.']
include_dirs = ['/home/mdgschmi/khinsen-mmtk-ae5b00d24d94/Include','.']
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

setup (name = "NoPotForceField",
       version = "1.0",
       description = "No forcefield term for MMTK",

       py_modules = ['NoPotFF'],
       ext_modules = [Extension('MMTK_NoPot',
                                ['MMTK_NoPot.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs,
                                libraries=libraries,
                                library_dirs=library_dirs)]
       )
