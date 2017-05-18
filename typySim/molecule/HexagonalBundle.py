from Molecule import Molecule
from ..geometry import shapes
import numpy as np
from math import ceil

class HexagonalBundle(Molecule):
  def __init__(self):
    super(HexagonalBundle,self).__init__() #Need to call parent class' constructor
    self.name = 'HexagonalBundle'
  # def build(self,nstems,nlayers,diameter=1.0,topType=0,bottomType=2,middleType=1):
  def build(self,bundle_radius,bundle_length,bead_diameter=1.0,topType=0,bottomType=2,middleType=1,base_gridsize=500):

    # pos = shapes.hexagonalspiral.positions(nstems,diameter)


    kwargs = {}
    kwargs['nx'] = base_gridsize
    kwargs['ny'] = base_gridsize
    kwargs['nz'] = 1
    kwargs['diameter'] = bead_diameter
    kwargs['topType']    = topType
    kwargs['middleType'] = topType
    kwargs['bottomType'] = topType
    kwargs['alternate_z'] = True
    molData = shapes.hexagonal.surface(**kwargs)

    pos = np.array([molData['x'], molData['y'], molData['z']]).T
    idex = pos.shape[0]/2
    idex = base_gridsize*base_gridsize/2 + base_gridsize/2
    cpos = pos[idex]
    dist = np.sqrt(np.sum(np.square(cpos - pos),axis=1))
    
    mask = ( dist<=bundle_radius )
    pos = pos[mask]
    pos -= cpos


    new_pos = []
    bonds = []
    types = [topType for i in range(pos.shape[0])]
    beadj = pos.shape[0]
    bundle_length = int(bundle_length)
    for pi,p in enumerate(pos):
      beadi = pi
      for l in range(1,bundle_length):
        new_pos.append([p[0],p[1],p[2]+l*bead_diameter])
        bonds.append([beadi,beadj])
        types.append(middleType)
        beadi = beadj
        beadj = beadj + 1
      types[-1] = bottomType

    pos = np.append(pos,new_pos,axis=0)
    types = np.array(types)

    

    
    molData = {}
    molData['x'] = pos[:,0]
    molData['y'] = pos[:,1]
    molData['z'] = pos[:,2]
    molData['types'] = types
    molData['bonds'] = bonds
    return molData,None
