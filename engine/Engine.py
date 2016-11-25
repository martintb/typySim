
class Engine(object):
  '''
  Base class for classes which procedurally update a :class:`System`
  positions, types, molecules, etc.

  .. Note:
    All subclasses should call super(`Class`,self).__init__() in their
    implementations of __init__(self).

  '''
  def __init__(self):
    self.system = None
    self.name = 'BaseEngine'
  def run(self):
    raise NotImplementedError('{} subclass of Engine has not implemented a run method!'.format(self.name))
