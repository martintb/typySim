from setuptools import setup,Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np

import typySim.Version
version = typySim.Version.short_version

'''
To compile the cython plugin, run this command:
python compile.py build_ext --inplace
'''
ext_modules = [
                Extension('*', 
                        [ 'typySim/core/cy/*.pyx' ],
												include_dirs=[np.get_include()],
 												extra_compile_args=['-fopenmp'],
 		     								extra_link_args=['-fopenmp']),
                Extension('*', 
                        [ 'typySim/compute/cy/*.pyx' ],
												include_dirs=[np.get_include()],
 												extra_compile_args=['-fopenmp','-Wno-maybe-uninitialized'],
 		     								extra_link_args=['-fopenmp','-Wno-maybe-uninitialized']),
                Extension('*', 
                        [ 'typySim/potential/cy/*.pyx' ],
												include_dirs=[np.get_include()],
 												extra_compile_args=['-fopenmp'],
 		     								extra_link_args=['-fopenmp']),
						  ]

setup(
    name='typySim',
	  ext_modules= cythonize(ext_modules,language='c++'),
    description='Python based molecular simulation engine',
    author='Tyler B. Martin',
    author_email = 'tyler.martin@nist.gov',
    version=version,
    packages=['typySim'],
    # scripts=['bin/typyViewer'],
    license='LICENSE.txt',
    long_description=open('README.md').read(),
)
