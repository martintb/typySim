#!python
# distutils: language=c++
# #cython: boundscheck=False
# #cython: wraparound=False
# #cython: cdivision=True
# #cython: nonecheck=False

import numpy as np
cimport numpy as np

from Compute cimport Compute
from libc.math cimport fabs as c_fabs
from libc.math cimport sqrt as c_sqrt
from libc.math cimport pow as c_pow
from libcpp.vector cimport vector

ctypedef double (*PotentialPtr)(double,double,double,double)


cdef double LennardJones(double r, double epsilon,double sigma, double rcut):
  cdef double U = 0
  cdef double r2inv,r6inv,r12inv,sig6,sig12;
  if r<rcut:
    r2inv = 1/(r*r)
    r6inv = r2inv*r2inv*r2inv
    r12inv = r6inv*r6inv
    sig6 = c_pow(sigma,6)
    sig12 = c_pow(sigma,12)
    U += 4*epsilon * (sig12*r12inv - sig6*r6inv)
  return U

cdef double HardSphere(double r, double epsilon,double sigma, double rcut):
  cdef double BIG = 1e9
  cdef double U = 0
  cdef double r2inv,r6inv,r12inv,sig6,sig12;
  if r<sigma:
    U = BIG
  return U
    
    

cdef class PotentialEnergy(Compute):
  cdef vector[ vector[PotentialPtr] ] PotentialMatrix;
  cdef double[:,:] epsilon_matrix 
  cdef double[:,:] sigma_matrix  
  cdef double[:,:] rcut_matrix  
  def __init__(self,system):
    super(PotentialEnergy,self).__init__()
    self.system = system
    self.build_matrices()
  def build_matrices(self):
    self.epsilon_matrix = np.array(self.system.PairTable.get_matrix('epsilon'))
    self.sigma_matrix   = np.array(self.system.PairTable.get_matrix('sigma'))
    self.rcut_matrix    = np.array(self.system.PairTable.get_matrix('rcut'))
    cdef long N = self.epsilon_matrix.shape[0]
    cdef Py_ssize_t i,j
    cdef vector[PotentialPtr]  temp;
    for i in range(N):
      temp.clear()
      for j in range(N):
        if (self.PairTable('potential',i,j) == 'HS'):
          temp.push_back(HardSphere)
        elif (self.PairTable('potential',i,j) == 'LJ'):
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
    cdef PotentialPtr UFunk

    for i in range(N-1):
      for j in range(i,N):

        dx = c_fabs(x[j] - x[i])
        dy = c_fabs(y[j] - y[i])
        dz = c_fabs(z[j] - z[i])
        dist = c_sqrt(dx*dx + dy*dy + dz*dz)

        ti = types[i]
        tj = types[i]
        rcut  = self.rcut_matrix[ti,tj]
        eps   = self.epsilon_matrix[ti,tj]
        sig   = self.sigma_matrix[ti,tj]
        UFunk = self.PotentialMatrix[ti][tj]
        U += UFunk(dist,eps,sig,rcut)

    return U
     
