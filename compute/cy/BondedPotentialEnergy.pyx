#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

from cython.parallel import parallel,prange

import numpy as np
cimport numpy as np

from libc.stdio cimport printf
from libc.math cimport sqrt as c_sqrt
from libc.math cimport pow as c_pow
from libcpp.vector cimport vector

from Compute cimport Compute

from typySim.core.cy.Box cimport * 
from typySim.core.cy.cymath cimport binary_search
from typySim.potential.cy.AllPotentials cimport *

cdef class BondedPotentialEnergy(Compute):
  def __init__(self,system):
    super(BondedPotentialEnergy,self).__init__()
    self._name = "BondedPotentialEnergy"
    self.system = system
    self.box = system.box
    self.build_matrices()
  def build_matrices(self):
    self.k_matrix = np.array(self.system.BondedTable.get_matrix('k'))
    self.r0_matrix   = np.array(self.system.BondedTable.get_matrix('r0'))
    cdef long N = self.k_matrix.shape[0]
    cdef Py_ssize_t i,j
    cdef vector[PotentialPointer]  temp;
    for i in range(N):
      temp.clear()
      for j in range(N):
        if (self.system.BondedTable['potential',i,j] == 'Harmonic'):
          temp.push_back(Harmonic)
        elif (self.system.BondedTable['potential',i,j] == 'FENE'):
          temp.push_back(FENE)
        else:
          raise ValueError('Bond type not recognized!')
      self.PotentialMatrix.push_back(temp)
  def compute(self,partial_indices=None,trial_move=False,**kwargs):
    cdef double U = -1.2345
    cdef double[:] x,trial_x
    cdef double[:] y,trial_y
    cdef double[:] z,trial_z
    cdef long[:] types,trial_types
    cdef long[:] indices
    cdef long[:,:] bonds,trial_bond_pairlist

    #regardless we need all of the system info
    x     = self.system.x
    y     = self.system.y
    z     = self.system.z
    types = self.system.types

    if trial_move:
      trial_x = self.system.trial_x
      trial_y = self.system.trial_y
      trial_z = self.system.trial_z
      trial_types = self.system.trial_types
      trial_bond_pairlist = self.system.trial_bond_pairlist
      if len(trial_bond_pairlist)>0:
        U = self.compute_trial_move(
                                        x,y,z,types,
                                        trial_x,trial_y,trial_z,trial_types,trial_bond_pairlist
                                       )
      else:
        U = 0
    elif self.system.bonds.nbonds>0:
      bonds = self.system.bonds.bonds
      if partial_indices is not None:
        indices = np.sort(partial_indices)
        U = self.compute_partial_system(indices,x,y,z,types,bonds)
      else:
        U = self.compute_full_system(x,y,z,types,bonds)
    else:
      #No bonds so no bond energy!
      U = 0

    return U
  cdef double compute_full_system(self, double[:] x, double[:] y, double[:] z, long[:] types, long[:,:] bonds) nogil:
    cdef double U = 0
    cdef Py_ssize_t bead_i,bead_j,bond_j
    cdef long N = bonds.shape[0]
    cdef double dx,dy,dz,dist
    cdef double k,r0
    cdef long ti,tj

    for bead_i in prange(N,nogil=True,schedule='guided'):
      bond_j = 0
      bead_j = bonds[bead_i,bond_j]
      while bead_j != -1:
        if bead_j>bead_i:
          dx = x[bead_j] - x[bead_i]
          dy = y[bead_j] - y[bead_i]
          dz = z[bead_j] - z[bead_i]

          dx = self.box.wrap_dx(dx)
          dy = self.box.wrap_dy(dy)
          dz = self.box.wrap_dz(dz)

          dist = c_sqrt(dx*dx + dy*dy + dz*dz)

          ti = types[bead_i]
          tj = types[bead_j]
          k    = self.k_matrix[ti,tj]
          r0 = self.r0_matrix[ti,tj]
          U += self.PotentialMatrix[ti][tj](dist,k,r0,0)
        bond_j = bond_j + 1
        bead_j = bonds[bead_i,bond_j]
    return U
  cdef double compute_partial_system(self, long[:] indices, double[:] x, double[:] y, double[:] z, long[:] types, long[:,:] bonds) nogil:
    cdef double U = 0
    cdef Py_ssize_t bead_i,bead_j,i,bond_j
    cdef double dx,dy,dz,dist
    cdef double k,r0
    cdef long ti,tj
    cdef long num_indices = indices.shape[0]

    # for bead_i in range(num_indices):
    for i in prange(num_indices,nogil=True,schedule='guided'):
      bead_i = indices[i]
      bond_j = 0
      bead_j = bonds[bead_i,bond_j]
      while bead_j != -1:
        
        # This if statement is supposed to guard against double counting. 
        #   Case 1: If j is in indices, we need to guard against double counting. Hence the
        #           j>i.
        #   Case 2: If j is not in indices, we count no matter what.
        if (bead_j>bead_i) or (binary_search(bead_j,indices) == -1):
          dx = x[bead_j] - x[bead_i]
          dy = y[bead_j] - y[bead_i]
          dz = z[bead_j] - z[bead_i]

          dx = self.box.wrap_dx(dx)
          dy = self.box.wrap_dy(dy)
          dz = self.box.wrap_dz(dz)

          dist = c_sqrt(dx*dx + dy*dy + dz*dz)

          ti = types[bead_i]
          tj = types[bead_j]
          k    = self.k_matrix[ti,tj]
          r0 = self.r0_matrix[ti,tj]
          U += self.PotentialMatrix[ti][tj](dist,k,r0,0)
        bond_j = bond_j + 1
        bead_j = bonds[bead_i,bond_j]
    return U
  cdef double compute_trial_move(self, 
                                    double[:] x, double[:] y, double[:] z, long[:] types,
                                    double[:] trial_x, double[:] trial_y, double[:] trial_z, long[:] trial_types, 
                                    long[:,:] trial_bond_pairlist) nogil:
    cdef double U = 0
    cdef Py_ssize_t bead_i,bead_j,bond_i
    cdef long N = trial_bond_pairlist.shape[0]
    cdef double dx,dy,dz,dist
    cdef double k,r0
    cdef long ti,tj
    cdef long nbeads_system = x.shape[0]
    cdef long nbeads_trial = trial_x.shape[0]
    cdef double x1,y1,z1
    cdef double x2,y2,z2

    for bond_i in prange(N,nogil=True,schedule='guided'):
      bead_i = trial_bond_pairlist[bond_i,0]
      bead_j = trial_bond_pairlist[bond_i,1]

      if bead_i>=nbeads_system:
        bead_i = bead_i - nbeads_system
        ti = trial_types[bead_i]
        x2= trial_x[bead_i]
        y2= trial_y[bead_i]
        z2= trial_z[bead_i]
      else:
        ti = types[bead_i]
        x2= x[bead_i]
        y2= y[bead_i]
        z2= z[bead_i]

      if bead_j>=nbeads_system:
        bead_j = bead_j - nbeads_system
        tj = trial_types[bead_j]
        x1= trial_x[bead_j]
        y1= trial_y[bead_j]
        z1= trial_z[bead_j]
      else:
        tj = types[bead_j]
        x1= x[bead_j]
        y1= y[bead_j]
        z1= z[bead_j]

      dx = x2 - x1
      dy = y2 - y1
      dz = z2 - z1

      dx = self.box.wrap_dx(dx)
      dy = self.box.wrap_dy(dy)
      dz = self.box.wrap_dz(dz)

      dist = c_sqrt(dx*dx + dy*dy + dz*dz)

      k    = self.k_matrix[ti,tj]
      r0 = self.r0_matrix[ti,tj]
      U += self.PotentialMatrix[ti][tj](dist,k,r0,0)
    return U
     
