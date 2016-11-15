from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np
'''
To compile the cython plugin, run this command:
python compile.py build_ext --inplace
'''
ext_modules = [
                Extension('*', ['typySim/cy/*.pyx'],
												include_dirs=[np.get_include()],
 												extra_compile_args=['-fopenmp'],
 		     								extra_link_args=['-fopenmp']),
												]
setup(
	cmdclass = {'build_ext': build_ext},
	ext_modules= cythonize(ext_modules),
  script_name='compile.py',
  script_args=['build_ext','--inplace']
)
