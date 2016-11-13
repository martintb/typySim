from Molecule import Molecule

class DummyMolecule(Molecule):
  def __init__(self):
    super(DummyMolecule,self).__init__()
    self.name = 'DummyMolecule'
  def build(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def check_indices(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def check_system(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def reset(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def update_indices(self,mapping):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def build(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def take_snapshot(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def restore_snapshot(self):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
  def append(self,indices):
    raise TypeError('This is placeholder molecule and should not be directly interacted with!')
