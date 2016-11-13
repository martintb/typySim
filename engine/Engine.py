
class Engine(object):
  def __init__(self):
    self.system = None
    self.name = 'BaseEngine'
  def run(self):
    raise NotImplementedError('{} subclass of Engine has not implemented a run method!'.format(self.name))
