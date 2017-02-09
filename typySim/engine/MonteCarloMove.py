
class MonteCarloMove(object):
  def __init__(self):
    self.engine = None
    self.system = None
    self.name   = 'BaseMonteCarloMove'
    self.call_count       = 0 
    self.attempt_count    = 0 
    self.accept_count     = 0 
    self.reset('')
  def __repr__(self):
    return self.__str__()
  def __str__(self):
    if self.attempt_count>0:
      rate1 = self.accept_count/float(self.attempt_count)
    else:
      rate1 = -1.0
    if self.call_count>0:
      rate2 = self.accept_count/float(self.call_count)
    else:
      rate2 = -1.0
    return '<{}/{:4.3f}/{:4.3f}>'.format(self.name,rate1,rate2)
  def reset(self,string_base):
    self.string           = string_base
    self.accept           = False
    self.technical_abort  = False
    self.Unew             = None
  def set_engine(self,engine):
    self.engine = engine
    self.system = engine.system
  def attempt(self):
    raise NotImplementedError('{} has not implemented \'attempt\' method!'.format(self.name))
  def attempt(self):
    self._attempt()
    self.call_count        += 1
    if not self.technical_abort:
      self.attempt_count   += 1
    if self.accept:
      self.accept_count    += 1
    


