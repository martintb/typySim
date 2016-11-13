from Molecule import Molecule
from ..geometry import shapes
import numpy as np
from math import ceil

class HexagonalSurface(Molecule):
  def __init__(self):
    super(HexagonalSurface,self).__init__() #Need to call parent class' constructor
    self.name = 'HexagonalSurface'
  def build(self,lx,ly,nz,diameter=1.0,topType=0,bottomType=2,middleType=1):
    #Find a close hexagonal grid for the requested box size
    nx,ny = shapes.hexagonal.position2Index(lx,ly)
    # nx and ny must be even for pbc to work
    nx = ceil(nx/2.0)*2.0 
    ny = ceil(ny/2.0)*2.0 
    # Find the actual box size based on the fitted nx and ny
    lx_fit,ly_fit= shapes.hexagonal.index2Position(nx,ny)
    
    print '==================================='
    print '>>> Box (lx,ly):'
    print '--> initial: ({},{})'.format(lx,ly)
    print '--> final:   ({},{})'.format(lx_fit,ly_fit)
    print '-----------------------------------'
    print '>>> Hex Grid (nx,ny):'
    print '--> size: ({},{})'.format(nx,ny)
    print '==================================='
    
    molData = shapes.hexagonal.surface(nx,ny,nz,diameter)
    return molData
