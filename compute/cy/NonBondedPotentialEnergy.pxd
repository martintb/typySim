#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

from libc.math cimport pow as c_pow
from libcpp.vector cimport vector

from Compute cimport Compute

from typySim.core.cy.Box cimport * 
from typySim.core.cy.CellList cimport *
from typySim.potential.cy.AllPotentials cimport *

cdef class NonBondedPotentialEnergy(Compute):
  cdef vector[ vector[PotentialPointer] ] PotentialMatrix;
  cdef double[:,:] epsilon_matrix 
  cdef double[:,:] sigma_matrix  
  cdef double[:,:] rcut_matrix  
  cdef Box box
  cdef CellList neighbor_list
  cdef double calc(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil
  cdef double calc_nlist(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil
