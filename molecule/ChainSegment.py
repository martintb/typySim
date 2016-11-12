from Molecule import Molecule
from ..geometry import shapes
import numpy as np

class ChainSegment(Molecule):
  def __init__(self):
    super(ChainSegment,self).__init__() #Need to call parent class' constructor
    self.name = 'ChainSegment'
  def build(self):
    molData = shapes.line(length,diameter)
    return molData
