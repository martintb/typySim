from Molecule import Molecule
from ..geometry import shapes
import numpy as np

import ipdb; ist = ipdb.set_trace

class ChainSegment(Molecule):
  def __init__(self):
    super(ChainSegment,self).__init__() #Need to call parent class' constructor
    self.name = 'ChainSegment'
    
    self.properties['topology']       = None
    self.properties['connected_to']   = {}
    self.properties['chain_ends']     = []
    self.properties['entangled_with'] = []
  def update_properties(self):
    chain_end1 =  self.indices[0]
    chain_end2 =  self.indices[-1]
    self.properties['chain_ends'] = [chain_end1,chain_end2]
    self.properties['connected_to'] = {}

    low_outer_bonds = self.get_outer_bonds(chain_end1,-1)
    if len(low_outer_bonds)==1:
      outer_index_low = low_outer_bonds[0][0]
      outer_mol_low = self.system.molecule_map[outer_index_low]
      self.properties['connected_to'][chain_end1] = {'index':outer_index_low,'molecule':outer_mol_low}
    elif len(low_outer_bonds)>1:
      raise ValueError('More than 1 `outer` bond found for ChainSegment')


    high_outer_bonds = self.get_outer_bonds(chain_end2,-1)
    if len(high_outer_bonds)==1:
      outer_index_high = high_outer_bonds[0][0]
      outer_mol_high = self.system.molecule_map[outer_index_high]
      self.properties['connected_to'][chain_end2] = {'index':outer_index_high,'molecule':outer_mol_high}
    elif len(high_outer_bonds)>1:
      raise ValueError('More than 1 `outer` bond found for ChainSegment')

    if len(self.properties['connected_to']) == 0:
      self.properties['topology'] = 'free'
      ist()
    elif len(self.properties['connected_to']) == 1:
      self.properties['topology'] = 'tail'
    elif len(self.properties['connected_to']) == 2 and (outer_mol_low is outer_mol_high):
      self.properties['topology'] = 'loop'
    elif len(self.properties['connected_to']) == 2 and (outer_mol_low is not outer_mol_high):
      self.properties['topology'] = 'tie'
    else:
      raise ValueError('ChainSegment has non_supported connectivity!')

  def get_outer_bonds(self,sys_index,new_index):
    bonds = []
    for bond_j in list(self.system.bonds.bonds[sys_index]):
      if bond_j==-1:
        break
      elif bond_j not in self.indices:
        bonds.append([bond_j,new_index])
    return bonds
  def build(self,**kwargs):
    molData = shapes.line(**kwargs)
    return molData
