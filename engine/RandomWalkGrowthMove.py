from MonteCarloMove import MonteCarloMove
import ..molecule
from geometry import linalg

import numpy as np
import random


class RandomWalkGrowthMove(MonteCarloMove):
  def __init__(self):
    super(RandomWalkGrowthMove,self).__init__() #must call parent class' constructor
    self.name='RandomWalkGrowthMove'
    self.grow_from_types = None
    self.end_type = None
  def attempt(self):
    # find bead to grow from
    types = self.system.types
    type_mask = np.in1d(types,self.grow_from_types)
    num_elements = type_mask.sum()
    choice = random.randint(0,num_elements)
    grow_from = self.system.positions[type_mask][choice]

    # generate random position
    newVec = 1.0*np.random.random((1,3)) - 0.5
    newPos = newVec - grow_from
    newPos = linalg.normalize(newVec - grow_from)

    # check for overlaps
    dr = self.system.positions - newPos
    dr = np.sqrt(np.sum(np.square(dr),axis=1))
    if np.any(dr<1.0):
      return False #reject!
    else:
      newMol = Bead()
      molData =  dict(positions=[newPos],types=2)
      self.system.add_molecule(newMol,build=False,**molData)
      return True


    # return True/False based on move acceptance


