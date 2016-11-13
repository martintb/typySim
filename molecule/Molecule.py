import numpy as np

class Molecule(object):
  def __init__(self):
    self.indices = None
    self.system = None
    self.positions = None
    self.types = None
    self.name = 'BaseMolecule'
    self.snapshot = None
  def check_indices(self):
    if self.indices is None:
      raise ValueError('Molecule named \"{}\" has no indices!'.format(self.name))
  def check_system(self):
    if self.system is None:
      raise ValueError('Molecule named \"{}\" has no system!'.format(self.name))
  def reset(self):
    self.check_system()
    self.check_indices()
    rawIndices = self.system.safe2raw[self.indices]
    mask = np.ones_like(self.system.positions,dtype=bool)
    mask[rawIndices] = False
    self.positions = np.ma.array(self.system.positions,mask=mask)
    self.positions._sharedmask = False
    self.types = np.ma.array(self.system.types,mask=mask[:,0])
    self.types._sharedmask = False
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
    if self.snapshot is None:
      raise ValueError('No snapshot to restore!')
    if (self.positions is None) or (self.indices is None) or (self.system is None):
      raise ValueError('Molecule named \"{}\" cannot restore snapshot without being fully specified!'.format(self.name))

    oldPos = self.snapshot['positions']
    newPos = self.positions
    if oldPos.mask != newPos.mask:
      raise ValueError('The mask has changed since the snapshot was taken. Cannot restore.')

    for od,nd,m in zip(oldPos,newPos,oldPos.mask):
      if not all(m):
        nd[:]=od
    self.snapshot = None