#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

from Compute cimport Compute
from libcpp.vector cimport vector
from typySim.core.cy.Box cimport * 
from typySim.potential.cy.AllPotentials cimport *

cdef class BondedPotentialEnergy(Compute):
  cdef vector[ vector[PotentialPointer] ] PotentialMatrix;
  cdef double[:,:] k_matrix 
  cdef double[:,:] r0_matrix  
  cdef Box box
  cdef double calc(self, double[:] x, double[:] y, double[:] z, long[:] types, long[:,:] bonds) nogil
     
