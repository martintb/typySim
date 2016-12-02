#!python
# distutils: language=c++
# #cython: boundscheck=False
# #cython: wraparound=False
# #cython: cdivision=True
# #cython: nonecheck=False

import numpy as np
cimport numpy as np

from libc.math cimport sqrt as c_sqrt
from libc.math cimport pow as c_pow
from libcpp.vector cimport vector

from Compute cimport Compute

from typySim.core.cy.Box cimport Box 
from typySim.core.cy.CellList cimport CellList
from typySim.potential.cy.AllPotentials cimport *


cdef class NonBondedPotentialEnergy(Compute):
  cdef vector[ vector[PotentialPointer] ] PotentialMatrix;
  cdef double[:,:] epsilon_matrix 
  cdef double[:,:] sigma_matrix  
  cdef double[:,:] rcut_matrix  
  cdef Box box
  cdef CellList neighbor_list
  def __init__(self,system):
    super(NonBondedPotentialEnergy,self).__init__()
    self._name = "NonBondedPotentialEnergy"
    self.system = system
    self.box = system.box
    self.neighbor_list = system.box.neighbor_list
    self.build_matrices()
  def build_matrices(self):
    self.epsilon_matrix = np.array(self.system.PairTable.get_matrix('epsilon'))
    self.sigma_matrix   = np.array(self.system.PairTable.get_matrix('sigma'))
    self.rcut_matrix    = np.array(self.system.PairTable.get_matrix('rcut'))
    cdef long N = self.epsilon_matrix.shape[0]
    cdef Py_ssize_t i,j
    cdef vector[PotentialPointer]  temp;
    for i in range(N):
      temp.clear()
      for j in range(N):
        if (self.system.PairTable['potential',i,j] == 'HS'):
          temp.push_back(HardSphere)
        elif (self.system.PairTable['potential',i,j] == 'LJ'):
          temp.push_back(LennardJones)
      self.PotentialMatrix.push_back(temp)
  def compute(self):
    cdef double U 
    cdef double[:] x 
    cdef double[:] y
    cdef double[:] z
    cdef long[:] types

    x     = self.system.x
    y     = self.system.y
    z     = self.system.z
    types = self.system.types

    U = self.calc_potential(x,y,z,types)
    return U

  cdef double calc_potential(self, double[:] x, double[:] y, double[:] z, long[:] types):
    cdef double U = 0
    cdef Py_ssize_t i,j
    cdef long N = x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj
    cdef PotentialPointer UFunk

    for i in range(N-1):
      for j in range(i+1,N):

        dx = x[j] - x[i]
        dy = y[j] - y[i]
        dz = z[j] - z[i]

        dx = self.box.wrap_dx(dx)
        dy = self.box.wrap_dy(dy)
        dz = self.box.wrap_dz(dz)

        dist = c_sqrt(dx*dx + dy*dy + dz*dz)

        ti = types[i]
        tj = types[j]
        rcut  = self.rcut_matrix[ti,tj]
        eps   = self.epsilon_matrix[ti,tj]
        sig   = self.sigma_matrix[ti,tj]
        UFunk = self.PotentialMatrix[ti][tj]
        U += UFunk(dist,eps,sig,rcut)
    return U
     
