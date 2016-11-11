from Molecule import Molecule
from ..geometry import shapes
import numpy as np

class HexagonalPlane(Molecule):
  def build(self,nx,ny,nz,diameter):
    molData = shapes.hexagonal_surface(nx,ny,nz,diameter,gen_types=True)
    return molData
