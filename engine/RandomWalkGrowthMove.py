from MonteCarloMove import MonteCarloMove
from .. import molecule
from ..core import Selection
from ..geometry import linalg
import pdb; st = pdb.set_trace

import numpy as np
import random


class RandomWalkGrowthMove(MonteCarloMove):
  def __init__(self):
    super(RandomWalkGrowthMove,self).__init__() #must call parent class' constructor
    self.name='RandomWalkGrowthMove'
    self.grow_from_types = None
    self.chain_end_type = None
    self.chain_middle_type = None
  def attempt(self):
    type_selection = Selection(self.system.types,self.grow_from_types)
    growth_index = type_selection.random_choice()
    growth_position = self.system.positions[growth_index]
    growth_type = self.system.types[growth_index]
    
    # generate random position
    newVec = linalg.normalize(np.random.random((1,3)) - 0.5)
    new_position = newVec + growth_position

    # check for overlaps
    dr = self.system.positions - new_position
    dist = self.system.box.wrap_distance(dr)
    if np.any(dist<1.0):
      return False #reject!
    else:
      print '======================================='
      print 'Successful Move!'
      print 'growth_position:',growth_position
      print 'growth_type:',growth_type
      print 'new_position',new_position
      print '======================================='
      molData = {}
      molData['positions'] = self.system.box.wrap_position(new_position)
      molData['types'] = self.chain_end_type

      #find if chain segment with the growth bead exists
      if (growth_type == self.chain_end_type):
        self.system.add_atoms(**molData)
        self.system.types[growth_index] = self.chain_middle_type
        ChainMol = self.system.molecule_map[growth_index]
        ChainMol.attach(self.system.natoms-1)
      else:
        ChainMol = molecule.ChainSegment()
        self.system.add_molecule(ChainMol,molData=molData)
      self.system.bonds[-1].add(growth_index)
      self.system.bonds[growth_index].add(self.system.natoms-1)

      return True




