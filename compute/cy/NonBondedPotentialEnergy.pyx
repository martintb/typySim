#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

from cython.parallel import parallel,prange

import numpy as np
cimport numpy as np

from libc.math cimport sqrt as c_sqrt
from libc.math cimport pow as c_pow
from libcpp.vector cimport vector

from Compute cimport Compute

from typySim.core.cy.Box cimport * 
from typySim.core.cy.CellList cimport *
from typySim.potential.cy.AllPotentials cimport *

cdef class NonBondedPotentialEnergy(Compute):
  # cdef vector[ vector[PotentialPointer] ] PotentialMatrix;
  # cdef double[:,:] epsilon_matrix 
  # cdef double[:,:] sigma_matrix  
  # cdef double[:,:] rcut_matrix  
  # cdef Box box
  # cdef CellList neighbor_list
  def __init__(self,system):
    super(NonBondedPotentialEnergy,self).__init__()
    self._name = "NonBondedPotentialEnergy"
    self.system = system
    self.box = system.box
    self.neighbor_list = system.box.neighbor_list
    self.build_matrices()
  def build_matrices(self):
    self.epsilon_matrix = np.array(self.system.NonBondedTable.get_matrix('epsilon'))
    self.sigma_matrix   = np.array(self.system.NonBondedTable.get_matrix('sigma'))
    self.rcut_matrix    = np.array(self.system.NonBondedTable.get_matrix('rcut'))
    cdef long N = self.epsilon_matrix.shape[0]
    cdef Py_ssize_t i,j
    cdef vector[PotentialPointer]  temp;
    for i in range(N):
      temp.clear()
      for j in range(N):
        if (self.system.NonBondedTable['potential',i,j] == 'HardSphere'):
          temp.push_back(HardSphere)
        elif (self.system.NonBondedTable['potential',i,j] == 'LennardJones'):
          temp.push_back(LennardJones)
        else:
          raise ValueError('Potential type not recognized!')
      self.PotentialMatrix.push_back(temp)
  def compute(self,ignore_neighbor_list=False):
    cdef double U = -1.2345
    cdef double[:] x 
    cdef double[:] y
    cdef double[:] z
    cdef long[:] types

    x     = self.system.x
    y     = self.system.y
    z     = self.system.z
    types = self.system.types

    if (not ignore_neighbor_list) and (self.box.neighbor_list is not None):
      if not self.box.neighbor_list.ready:
        raise ValueError('The neighbor list is reporting that it is not ready!')
      U = self.calc_nlist(x,y,z,types)
    else:
      U = self.calc(x,y,z,types)
    self.values.append(U)
    return U
  cdef double calc(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil:
    cdef double U = 0
    cdef Py_ssize_t i,j
    cdef long N = x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj

    for i in prange(N-1,nogil=True,schedule='guided'):
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
        rcut    = self.rcut_matrix[ti,tj]
        epsilon = self.epsilon_matrix[ti,tj]
        sigma   = self.sigma_matrix[ti,tj]
        U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)
    return U
  cdef double calc_nlist(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil:
    cdef double U = 0
    cdef Py_ssize_t i,j
    cdef long N = x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj
    cdef long cellNo,currCell,cellNeighNo

    for i in prange(N-1,nogil=True,schedule='guided'):
      cellNo = self.neighbor_list.bead_cells[i] #get the cell number of the current bead
      for cellNeighNo in range(27): #loop over the "neighbor cells" of this cell
        currCell = self.neighbor_list.cell_neighs[cellNo,cellNeighNo] 
        j = self.neighbor_list.top[currCell] #get the highest indexed bead in the neighbor chain of this cell
        while j!=-1:
          if j>i: #We don't want to double count
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            dz = z[j] - z[i]

            dx = self.box.wrap_dx(dx)
            dy = self.box.wrap_dy(dy)
            dz = self.box.wrap_dz(dz)

            dist = c_sqrt(dx*dx + dy*dy + dz*dz)

            ti = types[i]
            tj = types[j]
            rcut    = self.rcut_matrix[ti,tj]
            epsilon = self.epsilon_matrix[ti,tj]
            sigma   = self.sigma_matrix[ti,tj]
            U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)
          j = self.neighbor_list.neigh[j]
    return U
     
