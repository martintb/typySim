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
  def compute(self):
    cdef double U = -1.2345
    cdef double[:] x 
    cdef double[:] y
    cdef double[:] z
    cdef long[:] types
    cdef long[:,:] bonds

    if self.system.bonds.nbonds>0:
      x     = self.system.x
      y     = self.system.y
      z     = self.system.z
      types = self.system.types
      bonds = self.system.bonds.bonds

      U = self.calc(x,y,z,types,bonds)
    else: 
      U = 0
    self.values.append(U)
    return U
  cdef double calc(self, double[:] x, double[:] y, double[:] z, long[:] types, long[:,:] bonds) nogil:
    cdef double U = 0
    cdef Py_ssize_t i,j,bond_j
    cdef long N = bonds.shape[0]
    cdef double dx,dy,dz,dist
    cdef double k,r0
    cdef long ti,tj

    # for i in prange(N-1,nogil=True,schedule='guided'):
    for i in range(N):
      bond_j = 0
      j = bonds[i,bond_j]
      while j != -1:
        if j<=i:
          dx = x[j] - x[i]
          dy = y[j] - y[i]
          dz = z[j] - z[i]

          dx = self.box.wrap_dx(dx)
          dy = self.box.wrap_dy(dy)
          dz = self.box.wrap_dz(dz)

          dist = c_sqrt(dx*dx + dy*dy + dz*dz)

          ti = types[i]
          tj = types[j]
          k    = self.k_matrix[ti,tj]
          r0 = self.r0_matrix[ti,tj]
          U += self.PotentialMatrix[ti][tj](dist,k,r0,0)
        bond_j = bond_j + 1
        j = bonds[i,bond_j]
    return U
     
