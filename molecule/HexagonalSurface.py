from Molecule import Molecule
from ..geometry import shapes
import numpy as np

class HexagonalSurface(Molecule):
  def __init__(self):
    super(HexagonalSurface,self).__init__() #Need to call parent class' constructor
    self.name = 'HexagonalSurface'
  def build(self,nx,ny,nz,diameter):
    molData = shapes.hexagonal_surface(nx,ny,nz,diameter,gen_types=True)
    return molData
