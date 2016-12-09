#!python
# distutils: language=c++

from Compute cimport Compute
from libcpp.vector cimport vector
from typySim.core.cy.Box cimport * 
from typySim.potential.cy.AllPotentials cimport *

cdef class BondedPotentialEnergy(Compute):
  cdef vector[ vector[PotentialPointer] ] PotentialMatrix;
  cdef double[:,:] k_matrix 
  cdef double[:,:] r0_matrix  
  cdef Box box
  cdef double compute_full_system(self, double[:] x, double[:] y, double[:] z, long[:] types, long[:,:] bonds) nogil
  cdef double compute_partial_system(self,long[:] indices, double[:] x, double[:] y, double[:] z, long[:] types, long[:,:] bonds) nogil
  cdef double compute_trial_move(self, 
                                 double[:] x, double[:] y, double[:] z, long[:] types,
                                 double[:] trial_x, double[:] trial_y, double[:] trial_z, long[:] trial_types, 
                                 long[:,:] trial_bond_pairlist) nogil
