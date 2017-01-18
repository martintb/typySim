from Molecule import Molecule
from ..geometry import shapes
import numpy as np

class ChainSegment(Molecule):
  def __init__(self):
    super(ChainSegment,self).__init__() #Need to call parent class' constructor
    self.name = 'ChainSegment'
    
    self.properties['topology']       = None
    self.properties['connected_to']   = {}
    self.properties['chain_ends']     = []
    self.properties['entangled_with'] = []
  def find_connected(self):
    pass
  def build(self,**kwargs):
    molData = shapes.line(**kwargs)
    return molData
