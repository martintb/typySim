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
from typySim.potential.cy.AllPotentials cimport *

from typySim.compute.cy.BondedPotentialEnergy cimport *
from typySim.compute.cy.NonBondedPotentialEnergy cimport *

cdef class TotalPotentialEnergy(Compute):
  cdef BondedPotentialEnergy BPE
  cdef NonBondedPotentialEnergy NBPE
  def __init__(self,system):
    super(TotalPotentialEnergy,self).__init__()
    self._name = "TotalPotentialEnergy"
    self.system = system
    self.BPE = BondedPotentialEnergy(system)
    self.NBPE = NonBondedPotentialEnergy(system)
  def compute(self,**kwargs):
    cdef double U_BPE = -1.2345
    cdef double U_NBPE = -1.2345


    U_BPE = self.BPE.compute(**kwargs)
    U_NBPE = self.NBPE.compute(**kwargs)


    return sum([U_NBPE,U_BPE])
     
