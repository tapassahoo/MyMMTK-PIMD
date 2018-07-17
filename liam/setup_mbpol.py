import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

#os.environ["CC"] = "g++" 
#os.environ["CXX"] = "g++"


# Here is how to use the library built above.
ext_modules=[
    Extension("mbpol_eval",
              sources = ["mbpol_eval.pyx"],
              language="c++",
	      include_dirs = [os.getcwd(),'/home/kpbishop/mbpol/mbpol-19-Dec-2013/mbpol2/',np.get_include()],  # path to .h file(s)
              library_dirs = [os.getcwd(),'/home/kpbishop/mbpol/mbpol-19-Dec-2013/mbpol2/'],  # path to .a or .so file(s)
#              extra_objects = ['../mbpol/libmbpol.a'],
	      libraries = ['/home/kpbishop/mbpol/mbpol-19-Dec-2013/mbpol2/mbpol'])
]

setup(
  name = 'Demos',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
)
