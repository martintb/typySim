from Molecule import Molecule
from ..geometry import shapes
from ..geometry import linalg
import numpy as np
from math import ceil
import ipdb; ist = ipdb.set_trace

class HexagonalChevron2(Molecule):
  def __init__(self):
    super(HexagonalChevron2,self).__init__() #Need to call parent class' constructor
    self.name = 'HexagonalChevron2'
  # def build(self,nstems,nlayers,diameter=1.0,topType=0,bottomType=2,middleType=1):
  def build(self,
            angle,
            nx,
            ny,
            nz,
            bead_diameter=1.0,
            topType=0,
            bottomType=2,
            middleType=1,
            base_gridsize=100):


    ########################
    ## BUILD BASE LATTICE ##
    ########################
    kwargs = {}
    kwargs['nx'] = nx
    kwargs['ny'] = ny
    kwargs['nz'] = nz
    kwargs['diameter'] = bead_diameter
    kwargs['topType']    = topType
    kwargs['middleType'] = middleType
    kwargs['bottomType'] = bottomType
    kwargs['alternate_z'] = True
    molData = shapes.hexagonal.surface(**kwargs)

    ########################
    ## ROTATE AND REFLECT ##
    ########################
    pos = np.array([molData[i] for i in ['x','y','z']]).T

    pivot_x = np.max(pos[:,0])
    pivot_y = np.max(pos[:,1])
    pivot_z = np.min(pos[:,2])
    pivot   = np.array([pivot_x,pivot_y,pivot_z])
    pos -= pivot
    bond_shift = pos.shape[0]

    Q1 = linalg.Quaternion(axis=[1.0,0.0,0.0],angle=angle)
    Q2 = linalg.Quaternion(qvector=[0.0,0.0,1.0,0.0])
    rot_pos1 = []
    rot_pos2 = []
    for p in pos:
      pnew1  = Q1.rotate(p)
      pnew2  = Q2.reflect(pnew1)
      rot_pos1.append(pnew1)
      rot_pos2.append(pnew2)

    #####################
    ## REMOVE OVERLAPS ##
    #####################
    rot_pos2 = np.add(rot_pos2,[0,bead_diameter,0])
    pos = np.append(rot_pos1,rot_pos2,axis=0) 
    types = np.append(molData['types'],molData['types'],axis=0) 
    bonds = np.array(molData['bonds'])
    bonds = np.append(bonds,bonds+bond_shift,axis=0) 


    molData = {}
    molData['x'] = pos[:,0]
    molData['y'] = pos[:,1]
    molData['z'] = pos[:,2]
    molData['types'] = types
    molData['bonds'] = bonds
    return molData,None
