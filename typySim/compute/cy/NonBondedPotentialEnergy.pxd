#!python
# distutils: language=c++

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
  cdef double compute_full_system_no_nlist(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil
  cdef double compute_full_system(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil
  cdef double compute_partial_system(self, long[:] indices, double[:] x, double[:] y, double[:] z, long[:] types) nogil
  cdef double compute_trial_move(self, 
                                 double[:] x, double[:] y, double[:] z, long[:] types,
                                 double[:] trial_x, double[:] trial_y, double[:] trial_z, long[:] trial_types) nogil
  cdef double compute_trial_move_with_ignored(self, long[:] indices,
                                 double[:] x, double[:] y, double[:] z, long[:] types,
                                 double[:] trial_x, double[:] trial_y, double[:] trial_z, long[:] trial_types) nogil
