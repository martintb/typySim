#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

from libc.math cimport cos as c_cos
from libc.math cimport pi as c_pi

cdef double Soft(double r, double epsilon, double sigma, double rcut) nogil:
  cdef double U = 0
  if r<rcut:
    U = epsilon * (1 + c_cos(c_pi*r/rcut))
  return U
