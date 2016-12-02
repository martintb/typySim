#!python
# distutils: language=c++
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False

cdef double LennardJones(double r, double epsilon,double sigma, double rcut)
