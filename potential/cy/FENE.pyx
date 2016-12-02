#!python
# distutils: language=c++
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False

from libc.math cimport log as c_log

cdef double FENE(double r, double k,double r0, double dummy) nogil:
  return 0.5*k*r0*r0*c_log(1-(r/r0)*(r/r0))
