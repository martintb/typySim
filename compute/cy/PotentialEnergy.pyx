#!python
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

ctypedef double (*uptr)(double,double,double,double)

cdef double LJ126(self,double r, double epsilon,double sigma, double rcut):
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

cdef double HS(self,double r, double epsilon,double sigma, double rcut):
  cdef double BIG = 1e9
  cdef double U = 0
  cdef double r2inv,r6inv,r12inv,sig6,sig12;
  if r<sigma:
    U = BIG
  return U
    
    

cdef class PotentialEnergy(Compute):
  cdef object PairTable
  def __init__(self,PairTable):
    super(PotentialEnergy,self).__init__()
  def compute(self):
    cdef double[:] x                = self.system.x
    cdef double[:] y                = self.system.y
    cdef double[:] z                = self.system.z
    cdef long[:] types              = self.system.types
    cdef double[:,:] epsilon_matrix = np.array(self.PairTable.get_matrix('epsilon'))
    cdef double[:,:] sigma_matrix   = np.array(self.PairTable.get_matrix('sigma'))
    cdef double[:,:] rcut_matrix    = np.array(self.PairTable.get_matrix('rcut'))
    cdef uptr T;
    cdef vector[uptr] Tt;


    self.calc_potential(x,y,z,types,epsilon_matrix,sigma_matrix,rcut_matrix)

  cdef double calc_potential(self,
                             double[:] x, 
                             double[:] y, 
                             double[:] z,
                             long[:] types,
                             double[:,:] epsilon_matrix,
                             double[:,:] sigma_matrix,
                             double[:,:] rcut_matrix):
    cdef double U = 0
    return U
     
