from Molecule import Molecule
from ..geometry import shapes
import numpy as np

class ChainSegment(Molecule):
  def __init__(self):
    super(ChainSegment,self).__init__() #Need to call parent class' constructor
    self.name = 'ChainSegment'
    self.properties['tail'] = False
    self.properties['loop'] = False
    self.properties['tie']  = False
  def build(self,**kwargs):
    molData = shapes.line(**kwargs)
    return molData
