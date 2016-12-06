import warnings
class DummyMolecule(object):
  def __init__(self):
    self.name      = 'DummyMolecule'
    self.indices   = None
    self.system    = None
    self.positions = None
    self.types     = None
    self.bonds     = None
    self.snapshot  = None
  def __repr__(self):
    return "<{}>".format(self.name)
  def __str__(self):
    return "<{}>".format(self.name)
  def build(self):
    raise warning.warn('==> molecule.build() was called on DummyMolecule...')
  def check_indices(self):
    raise warning.warn('==> molecule.check_indices() was called on DummyMolecule...')
  def check_system(self):
    raise warning.warn('==> molecule.check_system() was called on DummyMolecule...')
  def reset(self,mapping):
    raise warning.warn('==> molecule.reset() was called on DummyMolecule...')
  def build(self):
    raise warning.warn('==> molecule.build() was called on DummyMolecule...')
  def take_snapshot(self):
    raise warning.warn('==> molecule.take_snapshot() was called on DummyMolecule...')
  def restore_snapshot(self):
    raise warning.warn('==> molecule.restore_snapshot() was called on DummyMolecule...')
  def attach(self,indices):
    import pdb; pdb.set_trace()
    raise warning.warn('==> molecule.attach() was called on DummyMolecule...')
  def detach(self,indices):
    import pdb; pdb.set_trace()
    raise warning.warn('==> molecule.detach() was called on DummyMolecule...')
  def isDummy(self):
    '''Check identity against global DummyMolecule sentinel.'''
    return True
    # if self is (self.system.DummyMolecule):
    #   return True
    # else:
    #   return False
