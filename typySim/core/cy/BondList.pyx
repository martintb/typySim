#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
import numpy as np
cimport numpy as np
from scipy.sparse.lil import lil_matrix
from scipy.sparse.csgraph import connected_components

intType = np.long
floatType = np.float32
doubleType = np.float
ctypedef np.int_t cIntType
ctypedef np.float32_t cFloatType
ctypedef np.float_t cDoubleType


cdef class BondList:
  cdef public long[:,:] bonds
  cdef public long max_bonds_perbead
  cdef public nbonds
  cdef bint init
  def __init__(self):
    self.max_bonds_perbead = 5
    self.init=False
    self.nbonds = 0
  def __getitem__(self,long index):
    return self.bonds[index]
  def reset(self,long num):
    self.init=False
    self.expand(num)
  def expand(self,long num):
    if self.init is False:
      self.bonds    = np.full((num,self.max_bonds_perbead),-1,dtype=intType)
      self.init = True
    else:
      self.bonds    = np.append(self.bonds, np.full((num,self.max_bonds_perbead),-1,dtype=intType),axis=0)
    self.nbonds = self.bonds.shape[0]
  def shrink(self,list indices):
    self.bonds  = np.delete(self.bonds,indices,axis=0)
  def add(self,i,j,shiftVal):
    '''
    Parameters
    ----------
    i,j : int, *required*
        Bond indices to be added to the :class:`System`.

        Bond indices passed as :class:`int` are assumed to be relative to the group of beads 
        being added.This means the passing [[0,1]] will add a bond between the first and 
        second beads **of the beads currently being added** and not the first and second 
        beads of the system. This effect is achieved by shifting the passed bond indices 
        by self.nbeads.

        Bond indices passed as :class:`str` indicate that you wish to bond against an index
        that has already been added o the system. These bond indices are made positive, but 
        left unshifted. 

        .. Warning::
          No sanity checks are carried out on the bond indices. Make sure you don't attempt
          to bond an bead index that doesn't exist in the system!

        .. Warning::
          This approach is hacky and bad and alternative ideas would be accepted. The basic problem
          is that a passed molecule needs to be able to specify all new internal bonds and any bonds
          to existing beads simultaneously. Bonds between new and old beads need to be possible. 
    shift : bool, *optional*
        Value to use as the shift-value if the beads are passed as :class:`int`. The default value
        is :func:`self.nbeads`.
     '''
    cdef long bondi,bondj

    if isinstance(i,basestring):
      bondi = int(i)
    else:
      bondi = i+shiftVal

    if isinstance(j,basestring):
      bondj = int(j)
    else:
      bondj = j+shiftVal

    self.insert_bond(bondi,bondj)
    self.insert_bond(bondj,bondi)
  def remove(self,long i,long j,long shiftVal):
    cdef long bondi,bondj

    if isinstance(i,basestring):
      bondi = int(i)
    else:
      bondi = i+shiftVal

    if isinstance(j,basestring):
      bondj = int(j)
    else:
      bondj = j+shiftVal

    self.excise_bond(bondi,bondj)
    self.excise_bond(bondj,bondi)
  cdef void insert_bond(self,long bondi,long bondj):
    cdef long new_j = bondj
    cdef long j,bond_num
    for bond_num in range(self.max_bonds_perbead):
      j = self.bonds[bondi,bond_num]

      if j == -1:
        self.bonds[bondi,bond_num] = new_j
        break
      elif j == new_j:
        break # This bond already exists! Abort!
      elif j > new_j:
        self.bonds[bondi,bond_num] = new_j
        new_j = j

      if bond_num==(self.max_bonds_perbead-1):
        raise ValueError('Too many bonds to bead {} when adding bead {}: {}'.format(bondi,bondj,self.bonds[bondi]))
  cdef void excise_bond(self,long bondi,long bondj):
    cdef long new_j = bondj
    cdef long j,bond_num

    #first loop over the list and change this bond to -1
    for bond_num in range(self.max_bonds_perbead):
      if self.bonds[bondi, bond_num] == bondj:
        self.bonds[bondi,bond_num] = -1
        break

    # We want to put the -1's at the end of the array. This accomplishes that by
    # making the sort function think they are infinity (np.inf). The neat thing is 
    # that the underlying values are unchanged so we get the desired behavior. 
    cdef long[:] blist = np.array(sorted(self.bonds[bondi],key=lambda x: np.inf if (x==-1) else x))
    self.bonds[bondi][:] = blist 
  def remap(self,long[:] mapping):
    # Processing the bonds is a bit of a pain. We need to remove all
    # references to removed indices while also mapping from old to 
    # new index numbers.
    cdef long i,j,old_j,new_j
    cdef long this_j, next_j
    cdef long num_bonds = self.bonds.shape[0]
    cdef long[:] blist 
    cdef long inf = int(1e9)

    for i in range(num_bonds):
      for j in range(self.max_bonds_perbead):
        old_j = self.bonds[i,j]
        if old_j != -1:
          self.bonds[i,j] = mapping[old_j]
      #Move any new -1's to the end
      for j in range(self.max_bonds_perbead-1):
        this_j = self.bonds[i,j]
        next_j = self.bonds[i,(j+1)]
        if this_j == -1:
          self.bonds[i,j]  = next_j
          self.bonds[i,(j+1)]  = this_j



      # for j in range(self.max_bonds_perbead):
      #   if self.bonds[i,j] != -1:
      #     self.bonds[i,j] = mapping[self.bonds[i,j]]

      # # We want to put the -1 at the end of the array. This accomplishes that by
      # # making the sort function think they are infinity (np.inf). The neat thing is 
      # # that the underlying values are unchanged so we get the desired behavior. 
      # blist = np.array(sorted(self.bonds[i],key=lambda x: np.inf if (x==-1) else x))
      # self.bonds[i][:] = blist 
  def connected(self):
    cdef long num_bonds = self.bonds.shape[0]
    cdef object lil = lil_matrix((num_bonds,num_bonds))
    cdef long bondj,bondi,j
    for bondi in range(num_bonds):
      for j in range(self.max_bonds_perbead):
        bondj = self.bonds[bondi,j]
        if (bondj!=-1) and (bondi<bondj):
          lil[bondi,bondj] = 1.0
    return connected_components(lil,directed=False)
  def get_pairlist(self,index=None):
    cdef list pairlist = []
    cdef long beadi,beadj,bondj,
    if index is None:
      for beadi in range(self.nbonds):
        for bondj in range(self.max_bonds_perbead):
          beadj = self.bonds[beadi,bondj]
          if beadj!=-1 and beadj>beadi:
            pairlist.append([beadi,beadj])
    else:
      beadi = index
      for bondj in range(self.max_bonds_perbead):
        beadj = self.bonds[beadi,bondj]
        if beadj!=-1:
          pairlist.append([beadi,beadj])
    return pairlist


