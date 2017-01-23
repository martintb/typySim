import numpy as np
import logging
from typySim.core.cy.cyutil import n_closest

from typySim.molecule.MoleculeCalculator import MoleculeCalculator
from typySim.molecule.MoleculeDistributor import MoleculeDistributor

class Molecule(object):
  '''
  Primary parent class for all reference molecule objects.

  In typySim a molecule is a "reference object" that references data which is actually 
  contained in a :class:`System` object. This effect is achieved via the 
  :class:`numpy.masked_array` module which allows masked references to be created to 
  existing arrays e.g. system.x, system.y, system.types,etc. While the implementation 
  of masked_array leaves much to be desired, it's as close as we're going to get to 
  pointer-based referencing ala C/C++. 

  The goal of this class is to provide a lightweight way to refer to a subset of beads
  in a :class:`System` while simultaneously providing easy and safe convenience functions
  to manipulate the subset of beads. 

  .. Note:
    All subclasses should call super(`Class`,self).__init__() in their
    implementations of __init__(self). 


  Attributes
  ----------
  indices : int list
      List of indices in the system object that this molecule references. Note
      that the the underlying data structure is a set and the list interface is
      served via a python property construct with defined setter. 

  system : object
      Reference to the system this molecule is contained in

  x,y,z,types : arrays
      Masked array references to arrays contained in system

  bonds : list of lists
      Currently not a reference but a copy of system.bonds. This will likely change

  name : str
      Each subclassed molecule with override this 

  '''
  def __init__(self):
    self._indices = None
    self.system = None
    self.x = None
    self.y = None
    self.z = None
    self.types = None
    self.bonds = None
    self.name = 'BaseMolecule'

    self.calculate = MoleculeCalculator(self)
    self.distribute = MoleculeDistributor(self)

    self.properties={}
    self.properties['center_of_mass'] = None
  @property
  def size(self):
    return len(self._indices)
  @property
  def indices(self):
    '''List interface for the molecules indices'''
    return self._indices
    # return list(self._indices) # set
  @indices.setter
  def indices(self,values):
    '''Parameter values must be iteratble.'''
    # self._indices = set(values) # set
    self._indices = list(values)
  def isDummy(self):
    '''Check identity against global DummyMolecule sentinel.'''
    return False
  def __repr__(self):
    return "<{}.{}>".format(self.name,self.size)
  def __str__(self):
    return "<{}.{}>".format(self.name,self.size)
  def check_indices(self):
    if self._indices is None:
      raise ValueError('Molecule named \"{}\" has no indices!'.format(self.name))
  def check_system(self):
    if self.system is None:
      raise ValueError('Molecule named \"{}\" has no system!'.format(self.name))
  def reset(self,index_mapping=None):
    '''(Re)build this molecule's data references to the system containers'''
    self.check_system()
    self.check_indices()

    # If atoms were added or removed from the sytem , we need to correctly 
    # map the old indices to their new values so that this molecule points 
    # to the correct bead
    if index_mapping is not None:
      new_indices = []
      for i in self.indices:
        if index_mapping[i] !=-1:
          new_indices.append(index_mapping[i])
      # self._indices = set(new_indices) # set
      self._indices = new_indices

    if not self._indices:
      # Molecule is empty! Signal for removal
      return False

    # Mask is FALSE where the atom belongs to this molecule
    mask = np.ones_like(self.system.x,dtype=bool)
    try:
      mask[self.indices] = False
    except IndexError:
      import ipdb; ipdb.set_trace()

    # These masked arrays are tricky, much care must be taken as
    # they don't exactly work the way you'd hope they would
    self.x = np.ma.array(self.system.x,mask=mask)
    self.y = np.ma.array(self.system.y,mask=mask)
    self.z = np.ma.array(self.system.z,mask=mask)
    self.imx = np.ma.array(self.system.imx,mask=mask)
    self.imy = np.ma.array(self.system.imy,mask=mask)
    self.imz = np.ma.array(self.system.imz,mask=mask)
    self.types = np.ma.array(self.system.types,mask=mask)

    # This is to suppress warning that have to do with creating submasks.
    # We should be safe with our usage. 
    self.x._sharedmask = False
    self.y._sharedmask = False
    self.z._sharedmask = False
    self.imx._sharedmask = False
    self.imy._sharedmask = False
    self.imz._sharedmask = False
    self.types._sharedmask = False

    self.bonds = []
    for idex in self.indices:
      self.bonds.append(self.system.bonds[idex])

    self.calculate.center_of_mass()

    return True
  def build(self):
    raise NotImplementedError('Molecule named \"{}\" hasn\'t defined a build method!'.format(self.name))
  def add_indices(self,indices):
    '''Add new indices to current molecule. Also updates the system.molecule_map'''
    for index in indices:
      # self._indices.add(index) # set
      self._indices.append(index)
      oldMolecule = self.system.molecule_map[index]
      self.system.molecule_map[index] = self
      # If oldMolecule is a not a placeholder, we need to remove 
      # references to this index from oldMolecule.
      if not oldMolecule.isDummy():
        oldMolecule.remove_indices([index])
    self.reset()
  def remove_indices(self,indices):
    ''' Removes indices from current molecule. Potentially updates the system.molecule_map'''
    for index in indices:
       #self._indices.remove(index) # set?
      self._indices.remove(index) 
      if self.system.molecule_map[index] is self:
        self.system.molecule_map[index] = self.system.DummyMolecule
    if not self.isDummy():
      self.reset()




