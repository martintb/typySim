from Molecule import Molecule
from ..geometry import shapes
import numpy as np
from math import ceil
import ipdb; ist = ipdb.set_trace

class HexagonalChevron1(Molecule):
  def __init__(self):
    super(HexagonalChevron1,self).__init__() #Need to call parent class' constructor
    self.name = 'HexagonalChevron1'
  # def build(self,nstems,nlayers,diameter=1.0,topType=0,bottomType=2,middleType=1):
  def build(self,
            arm_angle,
            arm_width,
            thickness,
            width,
            height,
            bead_diameter=1.0,
            topType=0,
            bottomType=2,
            middleType=1,
            base_gridsize=500):

    # pos = shapes.hexagonalspiral.positions(nstems,diameter)


    ########################
    ## BUILD BASE LATTICE ##
    ########################
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

    ###################
    ## TRACE CHEVRON ##
    ###################
    pos = np.array([molData['x'], molData['y'], molData['z']]).T
    idex1 = base_gridsize*base_gridsize/2 + base_gridsize/2
    cpos1 = pos[idex1]

  
    i,j,_ = shapes.hexagonal.position2Index(cpos1[0],cpos1[1]+arm_width)
    idex2 = i + j*base_gridsize + idex1
    cpos2 = pos[idex2]
    print cpos1,cpos2

    lineFunkLo = lambda x,th: np.tan(th)**(-1.0) * (x-cpos1[0]) + cpos1[1]
    lineFunkHi = lambda x,th: np.tan(th)**(-1.0) * (x-cpos2[0]) + cpos2[1]

    mask = np.ones_like(pos[:,0],dtype=bool)

    maskLeft  = pos[:,0] <= cpos1[0]
    mask1 = pos[maskLeft][:,1]>=lineFunkLo(pos[maskLeft][:,0],arm_angle)
    mask2 = pos[maskLeft][:,1]<=lineFunkHi(pos[maskLeft][:,0],arm_angle)
    mask[maskLeft]  = (mask1 & mask2)

    maskRight  = pos[:,0] >= cpos1[0]
    mask1 = pos[maskRight][:,1]>=lineFunkLo(pos[maskRight][:,0],-arm_angle)
    mask2 = pos[maskRight][:,1]<=lineFunkHi(pos[maskRight][:,0],-arm_angle)
    mask[maskRight]  = (mask1 & mask2)

    mask1 = abs(pos[:,0] - cpos1[0])<=width
    mask2 = abs(pos[:,1] - cpos1[1])<=height
    mask = mask & mask1 & mask2

    pos = pos[mask]
    pos -= cpos1

    ###############
    ## ADD DEPTH ##
    ###############
    new_pos = []
    bonds = []
    types = [topType for i in range(pos.shape[0])]
    beadj = pos.shape[0]
    thickness = int(thickness)
    for pi,p in enumerate(pos):
      beadi = pi
      for l in range(1,thickness):
        new_pos.append([p[0],p[1],p[2]+l*bead_diameter])
        bonds.append([beadi,beadj])
        types.append(middleType)
        beadi = beadj
        beadj = beadj + 1
      types[-1] = bottomType

    new_pos.append(cpos1)
    new_pos.append(cpos2)
    types.append(3)
    types.append(3)
    pos = np.append(pos,new_pos,axis=0)
    types = np.array(types)

    molData = {}
    molData['x'] = pos[:,0]
    molData['y'] = pos[:,1]
    molData['z'] = pos[:,2]
    molData['types'] = types
    molData['bonds'] = bonds
    return molData,None
