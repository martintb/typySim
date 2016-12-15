#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

cdef double HardSphere(double r, double epsilon,double sigma, double rcut) nogil:
  cdef double U = 0
  cdef double BIG = 1e4
  cdef double tol = 1e-6
  # binary, machine precision, and hexagonal lattices are dumb 
  if (sigma-r)>tol:
    U=BIG
  return U
