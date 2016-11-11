import numpy as np

class Molecule(object):
  def __init__(self):
    self.indices = None
    self.system = None
    self.positions = None
    self.name = 'BaseMolecule'
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
  def update_indices(self,mapping):
    raise NotImplementedError('This function needs to be implemented!')
  def build(self):
    raise NotImplementedError('Molecule named \"{}\" hasn\'t defined a build method!'.format(self.name))
