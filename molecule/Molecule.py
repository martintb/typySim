import numpy as np
import logging

class Molecule(object):
  def __init__(self):
    self.indices = None
    self.system = None
    self.x = None
    self.y = None
    self.z = None
    self.types = None
    self.bonds = None
    self.name = 'BaseMolecule'
    self.snapshot = None
    self.logger = logging.getLogger(__name__)
  def __repr__(self):
    return "<{}>".format(self.name)
  def __str__(self):
    return "<{}>".format(self.name)
  def check_indices(self):
    if self.indices is None:
      raise ValueError('Molecule named \"{}\" has no indices!'.format(self.name))
  def check_system(self):
    if self.system is None:
      raise ValueError('Molecule named \"{}\" has no system!'.format(self.name))
  def reset(self,index_mapping=None):
    self.check_system()
    self.check_indices()

    # If atoms were added or removed from the sytem , we need to correctly 
    # map the old indices to their new values so that this molecule points 
    # to the correct bead
    if index_mapping is not None:
      new_indices = []
      for i in self.indices:
        if index_mapping[i] is not None:
          new_indices.append(index_mapping[i])
      self.indices = set(new_indices)

    # Mask is FALSE where the atom belongs to this molecule
    mask = np.ones_like(self.system.x,dtype=bool)
    mask[list(self.indices)] = False

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
  def update_indices(self,mapping):
    raise NotImplementedError('This function needs to be implemented!')
  def build(self):
    raise NotImplementedError('Molecule named \"{}\" hasn\'t defined a build method!'.format(self.name))
  def take_snapshot(self):
    if (self.positions is None) or (self.indices is None) or (self.system is None):
      raise ValueError('Molecule named \"{}\" cannot take snapshot without being fully specified!'.format(self.name))
    self.snapshot = {}
    self.snapshot['positions'] = self.positions.copy()
  def restore_snapshot(self):
    raise NotImplementedError('This function needs to be implemented!')

  def attach(self,index):
    try:
      self.indices.add(index)
    except TypeError:
      raise TypeError('==> molecule.attach can only add one index at time!')
    oldMolecule = self.system.molecule_map[index]
    self.system.molecule_map[index] = self
    if oldMolecule is not self.system.DummyMolecule:
      oldMolecule.detach(index)
    self.reset()
  def detach(self,index):
    try:
      self.indices.remove(index)
    except TypeError:
      raise TypeError('==> molecule.detach can only remove one index at time!')
    if self.system.molecule_map[index] is self:
      self.system.molecule_map[index] = self.system.DummyMolecule
    if self is not self.system.DummyMolecule:
      self.reset()


