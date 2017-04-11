from math import sqrt
import numpy as np


def index2Position(i,j,k=0,r=0.5,alternate_z=False):
  '''
  i,j,k = x,y,z bead indices
  r = surface bead radius
  '''
  x = (2*i+((j)%2))*r
  y = sqrt(3)*(j)*r
  if alternate_z:
    z = 2.0 * k * r - (j%2)*r
  else:
    z = 2.0 * k * r
  return (x,y,z)

def position2Index(x,y,z=0,r=0.5,alternate_z=False):
  '''
  x,y = bead position
  r = surface bead radius
  '''
  j = int(float(y)/(sqrt(3)*r))
  i = int(0.5*(float(x)/r - ((j)%2)))
  if alternate_z:
    k = int(float(z-(j%2)*r)/r/2.0)
  else:
    k = int(float(z)/r/2.0)
  return (i,j,k)

def surface(nx,ny,nz,diameter,
            shift=True,
            topType=0,
            middleType=1,
            bottomType=2,
            alternate_z=False
            ):
  radius = diameter/2.0                    
  positions = []                           
  types = []
  chains = {}
  # bonds = []

  rowSize = nx
  layerSize = nx*ny
  n1d = nx
  n2d = nx*ny

  index3D = 0
  for k in range(nz):                      
    index2D = 0
    for j in range(ny):                    
      for i in range(nx):                  
        x,y,z = index2Position(i,j,k,r=radius,alternate_z=alternate_z)
        positions.append([x,y,z])

        if k==0:
          types.append(topType)
        elif k==(nz-1):
          types.append(bottomType)
        else:
          types.append(middleType)

        # if k>0:
        #   bead1 = i + j*n1d + (k-0)*n2d
        #   bead2 = i + j*n1d + (k-1)*n2d
        #   bonds.append([bead1,bead2])
        try:
          chains[index2D].append(index3D)
        except KeyError:
          chains[index2D] = [index3D]
        index2D += 1
        index3D += 1
  positions = np.array(positions)

  bonds = []
  for k,v in chains.items():
    bonds.extend([[i,j] for i,j in zip(v[:-1],v[1:])])

  if shift:
    positions -= np.average(positions,axis=0)
  molData ={}
  molData['x'] = positions[:,0]
  molData['y'] = positions[:,1]
  molData['z'] = positions[:,2]
  molData['types'] = types
  molData['bonds'] = bonds
  molData['chains'] = chains
  return molData
