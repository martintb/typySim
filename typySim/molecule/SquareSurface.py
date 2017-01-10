from Molecule import Molecule
from ..geometry import shapes
import numpy as np
from math import ceil

class SquareSurface(Molecule):
  def __init__(self):
    super(SquareSurface,self).__init__() #Need to call parent class' constructor
    self.name = 'SquareSurface'
  def build(self,lx,ly,nz,diameter=1.0,topType=0,bottomType=2,middleType=1):
    nx = int(lx/float(diameter))
    ny = int(ly/float(diameter))

    lx_fit = diameter*nx
    ly_fit = diameter*ny

    
    kwargs = {}
    kwargs['nx'] = nx
    kwargs['ny'] = ny
    kwargs['nz'] = nz
    kwargs['diameter'] = diameter
    kwargs['topType']    = topType
    kwargs['middleType'] = middleType
    kwargs['bottomType'] = bottomType
    molData = shapes.square.surface(**kwargs)
    boxData = {}
    boxData['lx'] = lx_fit
    boxData['ly'] = ly_fit
    boxData['lz'] = max(molData['x']) - min(molData['y'])+diameter
    return molData,boxData
