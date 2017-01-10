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
    nx = int(ceil(nx/2.0)*2.0)
    ny = int(ceil(ny/2.0)*2.0)
    # Find the actual box size based on the fitted nx and ny
    lx_fit,ly_fit= shapes.hexagonal.index2Position(nx,ny)
    

    logStr = ''' HexagonalSurface Created!
    Initial Box (lx,ly): ({},{}) 
    Final Box (lx,ly):   ({},{}) 
    Hex Grid (nx,ny):    ({},{})'''.format(lx,ly,lx_fit,ly_fit,nx,ny)
    # self.logger.debug(logStr)
    
    kwargs = {}
    kwargs['nx'] = nx
    kwargs['ny'] = ny
    kwargs['nz'] = nz
    kwargs['diameter'] = diameter
    kwargs['topType']    = topType
    kwargs['middleType'] = middleType
    kwargs['bottomType'] = bottomType
    molData = shapes.hexagonal.surface(**kwargs)
    boxData = {}
    boxData['lx'] = lx_fit
    boxData['ly'] = ly_fit
    boxData['lz'] = max(molData['x']) - min(molData['y'])
    return molData,boxData
