from random import choice
from numpy import in1d,arange,where

class Selection(object):
  def __init__(self,system):
    self.system = system
  def random_index(self,types=None):
    if types is not None:
      indices = where(in1d(self.system.types,types))[0]
    else:
      indices = arange(self.system.nbeads)
    return choice(indices)
  def random_molecule(self,names=[None]):
    molecules = []
    for name in names:
      if name is None:
        molecules.extend(self.system.molecules)
      else:
        molecules.extend(self.system.molecule_types[name])
    return choice(molecules)



