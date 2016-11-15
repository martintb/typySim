from MonteCarloMove import MonteCarloMove
from .. import molecule
from ..core import Selection
from ..geometry import linalg
import pdb; st = pdb.set_trace
import numpy as np
import random
import logging


class RandomWalkGrowthMove(MonteCarloMove):
  '''
  Grow non-overlapping random walk-based chains
  '''
  def __init__(self):
    super(RandomWalkGrowthMove,self).__init__() #must call parent class' constructor
    self.name='RandomWalkGrowthMove'
    self.grow_from_types = None
    self.chain_end_type = None
    self.chain_middle_type = None
    self.logger = logging.getLogger(__name__)
  def attempt(self):
    type_selection = Selection(self.system.types,self.grow_from_types)
    growth_index = type_selection.random_choice()
    growth_position = self.system.positions[growth_index]
    growth_type = self.system.types[growth_index]
    
    # generate random position
    newVec = linalg.normalize(np.random.random(3) - 0.5)
    new_position = newVec + growth_position

    # check for overlaps
    dr = self.system.positions - new_position
    dist = self.system.box.wrap_distance(dr)
    dist_mask = (dist<1.0)
    if np.any(dist_mask):
      return False #reject!
    else:
      self.logger.debug('''Accepted!
      growth_type:      {}
      growth_position:  {}
      new_position:     {}'''.format(growth_type,growth_position,new_position))

      molData = {}
      molData['positions'] = [self.system.box.wrap_position(new_position)]
      molData['types'] = self.chain_end_type
      new_index = self.system.nbeads

      if (growth_type == self.chain_end_type): # attaching bead to end of chain
        #add new bead to system
        self.system.add_beads(**molData) 

        #end of chain is no longer a growth type
        self.system.types[growth_index] = self.chain_middle_type

        #add new index to chain molecule
        self.system.molecule_map[growth_index].attach(self.system.nbeads-1)
      else: # attaching bead to surface

        #we need a new chain molecule to begin growing
        NewChainSegment = molecule.ChainSegment()

        #add new molecule to system
        self.system.add_molecule(NewChainSegment,molData=molData)

      #bonds between the surface/chain and new bead need to be added
      self.system.add_bond(new_index,growth_index)

      return True




