#!python
# distutils: language=c++

from typySim.core.cy.Box cimport * 

cdef long binary_search(long, long[:]) nogil
cpdef n_closest(long, double,double,double,double[:],double[:],double[:],Box)
