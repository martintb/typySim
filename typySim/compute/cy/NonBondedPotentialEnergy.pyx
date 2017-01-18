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
from typySim.core.cy.cyutil cimport binary_search

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
  def compute(self,partial_indices=None,trial_move=False,ignore_neighbor_list=False,ntrials=1,**kwargs):
    cdef double U = -1.2345
    cdef list Ulist = []
    cdef double[:] x,trial_x
    cdef double[:] y,trial_y
    cdef double[:] z,trial_z
    cdef long[:] types,trial_types
    cdef long[:] indices
    cdef long trial_num

    if (trial_move or (partial_indices is not None)) and ignore_neighbor_list:
      raise ValueError('Trial Move or Partial PE calculation is only supported with a neighbor_list!')

    x     = self.system.x
    y     = self.system.y
    z     = self.system.z
    types = self.system.types

    if trial_move:
      for trial_num in range(ntrials):
        trial_x = self.system.trial_x[trial_num]
        trial_y = self.system.trial_y[trial_num]
        trial_z = self.system.trial_z[trial_num]
        trial_types = self.system.trial_types[trial_num]
        if (partial_indices is not None):
          indices = np.sort(partial_indices)
          U = self.compute_trial_move_with_ignored(
                                                   indices,
                                                   x,y,z,types,
                                                   trial_x,trial_y,trial_z,trial_types
                                                  )
        else:
          U = self.compute_trial_move(
                                      x,y,z,types,
                                      trial_x,trial_y,trial_z,trial_types
                                     )
        Ulist.append(U)
    elif partial_indices is not None:
      indices = np.sort(partial_indices)
      Ulist.append(self.compute_partial_system(indices,x,y,z,types))
    else:
      if (not ignore_neighbor_list) and (self.box.neighbor_list is not None):
        if not self.box.neighbor_list.ready:
          raise ValueError('The neighbor list is reporting that it is not ready!')
        U = self.compute_full_system(x,y,z,types)
        Ulist.append(U)
      else:
        U = self.compute_full_system_no_nlist(x,y,z,types)
        Ulist.append(U)

    return Ulist
  cdef double compute_full_system_no_nlist(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil:
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
  cdef double compute_full_system(self, double[:] x, double[:] y, double[:] z, long[:] types) nogil:
    cdef double U = 0
    cdef Py_ssize_t i,j
    cdef long N = x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj
    cdef long cellNo,currCell,cellNeighNo

    for i in prange(N,nogil=True,schedule='guided'):
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
  cdef double compute_partial_system(self, long[:] indices, double[:] x, double[:] y, double[:] z, long[:] types) nogil:
    cdef double U = 0
    cdef Py_ssize_t i,j,bead_i,bead_j
    cdef long N = x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj
    cdef long cellNo,currCell,cellNeighNo
    cdef long num_indices = indices.shape[0]

    for i in prange(num_indices,nogil=True,schedule='guided'):
      bead_i = indices[i]
      cellNo = self.neighbor_list.bead_cells[bead_i] #get the cell number of the current bead
      for cellNeighNo in range(27): #loop over the "neighbor cells" of this cell
        currCell = self.neighbor_list.cell_neighs[cellNo,cellNeighNo] 
        bead_j = self.neighbor_list.top[currCell] #get the highest indexed bead in the neighbor chain of this cell
        while bead_j!=-1:
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
            rcut    = self.rcut_matrix[ti,tj]
            epsilon = self.epsilon_matrix[ti,tj]
            sigma   = self.sigma_matrix[ti,tj]
            U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)
          bead_j = self.neighbor_list.neigh[bead_j]
    return U
  cdef double compute_trial_move(self, 
                                 double[:] x, double[:] y, double[:] z, long[:] types,
                                 double[:] trial_x, double[:] trial_y, double[:] trial_z, long[:] trial_types) nogil:
    cdef double U = 0
    cdef Py_ssize_t bead_i,bead_j,i
    cdef long N_trial = trial_x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj
    cdef long cellNo,currCell,cellNeighNo
    cdef long ix,iy,iz
    cdef double x1,y1,z1
    cdef double x2,y2,z2
    cdef double ndx = self.neighbor_list.dx
    cdef double ndy = self.neighbor_list.dy
    cdef double ndz = self.neighbor_list.dz
    cdef double bx = self.neighbor_list.bx
    cdef double by = self.neighbor_list.by
    cdef double bz = self.neighbor_list.bz

    # for bead_i in range(N_trial-1):
    for bead_i in prange(N_trial,nogil=True,schedule='guided'):
      x1 = trial_x[bead_i]
      y1 = trial_y[bead_i]
      z1 = trial_z[bead_i]

      #intra
      for bead_j in range(bead_i+1,N_trial):

        dx = trial_x[bead_j] - x1
        dy = trial_y[bead_j] - y1
        dz = trial_z[bead_j] - z1

        dx = self.box.wrap_dx(dx)
        dy = self.box.wrap_dy(dy)
        dz = self.box.wrap_dz(dz)

        dist = c_sqrt(dx*dx + dy*dy + dz*dz)

        ti = trial_types[bead_i]
        tj = trial_types[bead_j]
        rcut    = self.rcut_matrix[ti,tj]
        epsilon = self.epsilon_matrix[ti,tj]
        sigma   = self.sigma_matrix[ti,tj]
        U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)

      #inter
      ix = self.neighbor_list.pos2idex(x1,ndx,bx,True)
      iy = self.neighbor_list.pos2idex(y1,ndy,by,True)
      iz = self.neighbor_list.pos2idex(z1,ndz,bz,True)
      cellNo = self.neighbor_list.idex2cell(ix,iy,iz)
      for cellNeighNo in range(27): #loop over the "neighbor cells" of this cell
        currCell = self.neighbor_list.cell_neighs[cellNo,cellNeighNo] 
        bead_j = self.neighbor_list.top[currCell] #get the highest indexed bead in the neighbor chain of this cell
        while bead_j!=-1:
          dx = x[bead_j] - x1
          dy = y[bead_j] - y1
          dz = z[bead_j] - z1

          dx = self.box.wrap_dx(dx)
          dy = self.box.wrap_dy(dy)
          dz = self.box.wrap_dz(dz)

          dist = c_sqrt(dx*dx + dy*dy + dz*dz)

          ti = trial_types[bead_i]
          tj = types[bead_j]
          rcut    = self.rcut_matrix[ti,tj]
          epsilon = self.epsilon_matrix[ti,tj]
          sigma   = self.sigma_matrix[ti,tj]
          U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)
          bead_j = self.neighbor_list.neigh[bead_j]
    return U
  cdef double compute_trial_move_with_ignored(self, long[:] ignored,
                                 double[:] x, double[:] y, double[:] z, long[:] types,
                                 double[:] trial_x, double[:] trial_y, double[:] trial_z, long[:] trial_types) nogil:
    cdef double U = 0
    cdef Py_ssize_t bead_i,bead_j,i
    cdef long N_trial = trial_x.shape[0]
    cdef double dx,dy,dz,dist
    cdef double epsilon,sigma,rcut
    cdef long ti,tj
    cdef long cellNo,currCell,cellNeighNo
    cdef long ix,iy,iz
    cdef double x1,y1,z1
    cdef double x2,y2,z2

    #intra
    # for bead_i in range(N_trial-1):
    for bead_i in prange(N_trial-1,nogil=True,schedule='guided'):
      for bead_j in range(bead_i+1,N_trial):

        dx = trial_x[bead_j] - trial_x[bead_i]
        dy = trial_y[bead_j] - trial_y[bead_i]
        dz = trial_z[bead_j] - trial_z[bead_i]

        dx = self.box.wrap_dx(dx)
        dy = self.box.wrap_dy(dy)
        dz = self.box.wrap_dz(dz)

        dist = c_sqrt(dx*dx + dy*dy + dz*dz)

        ti = trial_types[bead_i]
        tj = trial_types[bead_j]
        rcut    = self.rcut_matrix[ti,tj]
        epsilon = self.epsilon_matrix[ti,tj]
        sigma   = self.sigma_matrix[ti,tj]
        U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)

    #inter
    # for bead_i in range(N_trial):
    for bead_i in prange(N_trial,nogil=True,schedule='guided'):
      x1 = trial_x[bead_i]
      y1 = trial_y[bead_i]
      z1 = trial_z[bead_i]
      ix = self.neighbor_list.pos2idex(x1,self.neighbor_list.dx,self.neighbor_list.bx,True)
      iy = self.neighbor_list.pos2idex(y1,self.neighbor_list.dy,self.neighbor_list.by,True)
      iz = self.neighbor_list.pos2idex(z1,self.neighbor_list.dz,self.neighbor_list.bz,True)
      cellNo = self.neighbor_list.idex2cell(ix,iy,iz)
      for cellNeighNo in range(27): #loop over the "neighbor cells" of this cell
        currCell = self.neighbor_list.cell_neighs[cellNo,cellNeighNo] 
        bead_j = self.neighbor_list.top[currCell] #get the highest indexed bead in the neighbor chain of this cell
        while bead_j!=-1:
          # This if statement is supposed to guard against double counting. 
          #   Case 1: If j is in indices, we need to guard against double counting. Hence the
          #           j>i.
          #   Case 2: If j is not in indices, we count no matter what.
          if binary_search(bead_j,ignored) == -1:
            dx = x[bead_j] - x1
            dy = y[bead_j] - y1
            dz = z[bead_j] - z1

            dx = self.box.wrap_dx(dx)
            dy = self.box.wrap_dy(dy)
            dz = self.box.wrap_dz(dz)

            dist = c_sqrt(dx*dx + dy*dy + dz*dz)

            ti = trial_types[bead_i]
            tj = types[bead_j]
            rcut    = self.rcut_matrix[ti,tj]
            epsilon = self.epsilon_matrix[ti,tj]
            sigma   = self.sigma_matrix[ti,tj]
            U += self.PotentialMatrix[ti][tj](dist,epsilon,sigma,rcut)
          bead_j = self.neighbor_list.neigh[bead_j]
    return U
     
