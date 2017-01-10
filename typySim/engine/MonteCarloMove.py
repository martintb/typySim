
class MonteCarloMove(object):
  def __init__(self):
    self.engine = None
    self.system = None
    self.name = 'BaseMonteCarloMove'
  def attempt(self):
    raise NotImplementedError('{} has not implemented \'attempt\' method!'.format(self.name))
  @staticmethod
  def counter(fn):
    def _wrapper(*args,**kwargs):
      _wrapper.call_count += 1
      result =  fn(*args,**kwargs)
      if result is True:
        _wrapper.accept_count += 1
      return result
    _wrapper.call_count  = 0 
    _wrapper.accept_count = 0 
    _wrapper.__name__ = fn.__name__
    return _wrapper
    


