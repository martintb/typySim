import numpy as np
import logging
from typySim.core.cy.cyutil import n_closest

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
    self.properties={}
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
    # if self is (self.system.DummyMolecule):
    #   return True
    # else:
    #   return False
  def __repr__(self):
    return "<{}>".format(self.name)
  def __str__(self):
    return "<{}>".format(self.name)
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
    self.types = np.ma.array(self.system.types,mask=mask)

    # This is to suppress warning that have to do with creating submasks.
    # We should be safe with our usage. 
    self.x._sharedmask = False
    self.y._sharedmask = False
    self.z._sharedmask = False
    self.types._sharedmask = False

    self.bonds = []
    for idex in self.indices:
      self.bonds.append(self.system.bonds[idex])
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
  def distribute_by_index(self,index_groups):
    '''
    Break this molecule up into smaller molecules of identical type.

    Parameters
    ----------
    index_groups : list, *Required*
        Index group is a two-dimensional list where len(index_groups) is 
        the number of sub-molecules to create and each list item is a list
        of indices to use for each molecule.

        .. Warning:
          All indices must be accounted for in the new molecules. An exception
          will be raised if this is not true.
    
    '''
    #Flatten list of lists into a 1-D set
    index_check = set([i for group in index_groups for i in group])
    # if index_check!=self._indices: #set
    if index_check!=set(self._indices):
      raise ValueError('Mismatch between passed indices and molecule indices.')

    new_molecules = []
    for group in index_groups:
      new_molecule = self.__class__()
      new_molecule.indices = group
      self.system.add_molecule(new_molecule)
      new_molecules.append(new_molecule)
    self.system.remove_molecule(molecule=self,remove_beads=False)
    return new_molecules
  def distribute_by_type_proximity(self,typeList,box=None):
    ''' Distribute a single molecule into multiple molecules based on proximity to a bead type

    '''
    if len(typeList)<2:
      raise ValueError('Need at least two types to distribute between!')

    if box is None:
      try:
        box = self.system.box
      except AttributeError:
        raise ValueError('Molecule must be connected to system, or explicit box object must be passed!')


    x = self.x.compressed()
    y = self.y.compressed()
    z = self.z.compressed()
    types = self.types.compressed()
    indices = np.array(self._indices)


    # index_groups will be passed to self.distribute_by_index. It must be initialized with the
    # beads which are pre-grouped because they match a grouping type
    index_groups = []
    for i,typeVal in enumerate(typeList):
      index_groups.append(list(indices[types==typeVal]))

    # Need to combine all typeGroups. These beads are already grouped
    # as they define the groups themselves
    mask = (types==typeList[0])
    for typeVal in typeList[1:]:
      mask = np.logical_or(mask,types==typeVal)

    # Beads to be grouped
    x1Array  = x[~mask]
    y1Array  = y[~mask]
    z1Array  = z[~mask]
    indices1 = indices[~mask]

    # Beads which are pre-grouped because they match a proximity type
    x2Array = x[mask]
    y2Array = y[mask]
    z2Array = z[mask]
    types2  = types[mask]

    for i,(x1,y1,z1) in enumerate(zip(x1Array,y1Array,z1Array)):
      idex,dist = n_closest(1,x1,y1,z1,x2Array,y2Array,z2Array,box)
      for j,typeVal in enumerate(typeList):
        if types2[idex[0]] == typeVal:
          index_groups[j].append(indices1[i])
          break

    return self.distribute_by_index(index_groups)
  def center_of_mass(self,box=None):
    if box is None:
      try:
        box = self.system.box
      except AttributeError:
        raise ValueError('Molecule must be connected to system, or explicit box object must be passed!')

    mask = self.x.mask
    ux = self.system.x[~mask]
    uy = self.system.y[~mask]
    uz = self.system.z[~mask]

    imx = self.system.imx[~mask]
    imy = self.system.imy[~mask]
    imz = self.system.imz[~mask]

    bx = self.system.box.lx
    by = self.system.box.ly
    bz = self.system.box.lz

    ux = ux - imx*bx
    uy = uy - imy*by
    uz = uz - imz*bz

    return np.mean(ux),np.mean(uy),np.mean(uz)




