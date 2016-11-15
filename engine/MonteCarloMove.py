
class MonteCarloMove(object):
  def __init__(self):
    self.engine = None
    self.system = None
    self.name = 'BaseMonteCarloMove'
  def attempt(self):
    raise NotImplementedError('{} has not implemented \'attempt\' method!'.format(self.name))
