from random import choice
from numpy import in1d,arange,where

class Selection(object):
  def __init__(self,system):
    self.system = system
  @profile
  def random_index(self,types=None):
    if types is not None:
      indices = where(in1d(self.system.types,types))[0]
    else:
      indices = arange(self.system.nbeads)
    return choice(indices)



