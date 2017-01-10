#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

cdef double Harmonic(double r, double k,double r0, double dummy) nogil:
  return 0.5 * k * (r-r0)*(r-r0)
