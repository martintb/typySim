#!python
# distutils: language=c++
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False

from libc.math cimport pow as c_pow

cdef double LennardJones(double r, double epsilon,double sigma, double rcut):
  cdef double U = 0
  cdef double r2inv,r6inv,r12inv,sig6,sig12;
  if r<rcut:
    r2inv = 1/(r*r)
    r6inv = r2inv*r2inv*r2inv
    r12inv = r6inv*r6inv
    sig6 = c_pow(sigma,6)
    sig12 = c_pow(sigma,12)
    U += 4*epsilon * (sig12*r12inv - sig6*r6inv)
  return U
