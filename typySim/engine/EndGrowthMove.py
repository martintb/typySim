from MonteCarloMove import *
from typySim import molecule
from typySim.geometry import linalg
import numpy as np
import random

import ipdb as ipdb; ist = ipdb.set_trace


class EndGrowthMove(MonteCarloMove):
  '''
  Append or remove beads to the ends of growing chains
  '''
  def __init__(self,surface_neutral_type,surface_growth_type,chain_end_type,chain_middle_type):
    super(EndGrowthMove,self).__init__() #must call parent class' constructor
    self.name='EndGrowthMove'
    self.growth_types = [surface_growth_type,chain_end_type]
    self.chain_end_type = chain_end_type
    self.chain_middle_type = chain_middle_type
    self.surface_neutral_type = surface_neutral_type
  def _attempt(self):
    self.reset('EndGrowth')
    Uold = self.engine.TPE_list[-1]

    growth_index = self.system.select.random_index(types=self.growth_types)
    growth_x     = self.system.x[growth_index]
    growth_y     = self.system.y[growth_index]
    growth_z     = self.system.z[growth_index]
    growth_type  = self.system.types[growth_index]
    growth_mol   = self.system.molecule_map[growth_index]
    
    # generate random position. The magic at the end of the line transforms the size (3,) array
    # to a size (3,1) array which is important for the wrapping step.
    newVec = linalg.normalize(np.random.random(3) - 0.5)[np.newaxis].T
    newVec[0] += growth_x
    newVec[1] += growth_y
    newVec[2] += growth_z
    (new_x,new_y,new_z),(new_imx,new_imy,new_imz) = self.system.box.wrap_positions(x=newVec[0],y=newVec[1],z=newVec[2])
    new_index = self.system.nbeads
    new_types = self.chain_end_type
    new_bonds = [[growth_index,new_index]]
    self.system.set_trial_move(
                                x=[new_x],
                                y=[new_y],
                                z=[new_z],
                                imx=[new_imx],
                                imy=[new_imy],
                                imz=[new_imz],
                                types=[[new_types]],
                                bonds=new_bonds)

    Unew = Uold + sum(self.engine.TPE.compute(trial_move=True))
    if Unew<=Uold:
      self.accept=True
    else:
      prob = np.exp(-(Unew-Uold))
      ranf = np.random.rand()
      if ranf<prob:
        self.accept=True
      else:
        self.accept=False

    if self.accept:
      self.system.append_trial_move()
      new_index = [new_index]
      if (growth_type == self.chain_end_type): # attaching bead to end of chain
        #end of chain is no longer a growth type
        self.system.types[growth_index] = self.chain_middle_type

        #add new index to chain molecule
        # self.system.molecule_map[growth_index].add_indices(new_index)
        growth_mol.add_indices(new_index)

        growth_mol.properties['chain_ends'].append(new_index[0])
        if growth_mol.size>2: #Needs to be two because we just added a bead
          growth_mol.properties['chain_ends'].remove(growth_index)

      else: # attaching bead to surface

        #end of chain is no longer a growth type
        self.system.types[growth_index] = self.surface_neutral_type

        mol_index = np.where(self.system.molecules==growth_mol)[0]
        #we need a new chain molecule to begin growing
        NewChainSegment = molecule.ChainSegment()
        NewChainSegment.indices = new_index
        NewChainSegment.properties['topology']                   = 'tail'
        NewChainSegment.properties['connected_to'][new_index[0]] = {'index':growth_index,'molecule':mol_index}
        NewChainSegment.properties['chain_ends']                 = [new_index[0]]
        self.system.add_molecule(NewChainSegment)

      # XXX This should probably be called after every success of this move, but
      #     will likely cause performance loss. An alternative would be to be rigorous
      #     about calling reset() on each molecule that we work with. This is only
      #     a problem when the number of beads in the system changes, as this event
      #     causes the masked arrays to become stale (i.e. incorrectly sized)
      self.system.reset_all_molecules()

    self.Unew   = Unew
    return




