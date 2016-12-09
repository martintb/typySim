from MonteCarloMove import *
from typySim import molecule
from typySim.core import Selection
from typySim.geometry import linalg
import numpy as np
import random
import logging

# import ipdb as pdb; 
# import os; eexit = os._exit
# st = pdb.set_trace


class EndGrowthMove(MonteCarloMove):
  '''
  Append or remove beads to the ends of growing chains
  '''
  def __init__(self,surface_growth_type,chain_end_type,chain_middle_type):
    super(EndGrowthMove,self).__init__() #must call parent class' constructor
    self.name='EndGrowthMove'
    self.growth_types = [surface_growth_type,chain_end_type]
    self.chain_end_type = chain_end_type
    self.chain_middle_type = chain_middle_type
  @MonteCarloMove.counter
  def attempt(self):
    Uold = self.engine.TPE_list[-1]

    type_selection = Selection(self.system.types,self.growth_types)
    growth_index = type_selection.random_choice()
    growth_x = self.system.x[growth_index]
    growth_y = self.system.y[growth_index]
    growth_z = self.system.z[growth_index]
    growth_type = self.system.types[growth_index]
    
    # generate random position
    newVec = linalg.normalize(np.random.random(3) - 0.5)
    new_x = newVec[0] + growth_x
    new_y = newVec[1] + growth_y
    new_z = newVec[2] + growth_z
    new_x,new_y,new_z = self.system.box.numpy_wrap_position(x=new_x,y=new_y,z=new_z)
    new_index = self.system.nbeads

    molData = {}
    molData['x'] = [new_x]
    molData['y'] = [new_y]
    molData['z'] = [new_z]
    molData['types'] = [self.chain_end_type]

    # add new bead to system
    if (growth_type == self.chain_end_type): # attaching bead to end of chain

      #add new bead to system
      self.system.add_beads(**molData) 

      #end of chain is no longer a growth type
      self.system.types[growth_index] = self.chain_middle_type

      #add new index to chain molecule
      self.system.molecule_map[growth_index].add_indices([new_index])
    else: # attaching bead to surface
      #we need a new chain molecule to begin growing
      NewChainSegment = molecule.ChainSegment()

      #add new molecule to system
      self.system.add_molecule(NewChainSegment,**molData)



    Unew = self.engine.TPE.compute(ignore_neighbor_list=True)
    if Unew<Uold:
      accept=True
    else:
      prob = np.exp(-(Unew-Uold))
      ranf = np.random.rand()
      if ranf<prob:
        accept=True
      else:
        accept=False

    if not accept:
      #bonds between the surface/chain and new bead need to be added
      # self.system.bonds.remove(new_index,growth_index,0)

      #Ugh...we need to put the system back the way it was
      if (growth_type == self.chain_end_type): # bead was attached to end of chain
        mapping = self.system.remove_beads([new_index]) 
      else: # bead was attached to surface
        self.system.remove_molecule(molecule=NewChainSegment,remove_beads=True)

    return accept




