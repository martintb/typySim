#!python
# distutils: language=c++
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False

ctypedef double (*PotentialPointer)(double,double,double,double)
# cdef double LennardJones(double r, double epsilon,double sigma, double rcut)
# cdef double HardSphere(double r, double epsilon,double sigma, double rcut)
