from math import sqrt
import numpy as np

def surface(nx,ny,nz,diameter,
            shift=True,
            topType=0,
            middleType=1,
            bottomType=2,
            ):
  radius = diameter/2.0                    
  positions = []                           
  types = []
  chains = {}
  index3D = 0
  for k in range(nz):                      
    index2D = 0
    for j in range(ny):                    
      for i in range(nx):                  
        x = (2*i + ((j+k)%2))*radius
        y = sqrt(3)*(j+(k%2)/3.0)*radius
        z = 2.0/3.0 * sqrt(6) * k * radius
        positions.append([x,y,z])

        if k==0:
          types.append(topType)
        elif k==(nz-1):
          types.append(bottomType)
        else:
          types.append(middleType)

        try:
          chains[(index2D,k%2)].append(index3D)
        except KeyError:
          chains[(index2D,k%2)] = [index3D]
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

def index2Position(i,j,k=0,r=0.5):
  '''
  i,j,k = x,y,z bead indices
  r = surface bead radius
  '''
  x = (2*i+((j+k)%2))*r
  y = sqrt(3)*(j+(k%2)/3.0)*r
  return (x,y)

def position2Index(x,y,k=0.0,r=0.5):
  '''
  x,y = bead position
  r = surface bead radius
  '''
  j = int(float(y)/(sqrt(3)*r) - (k%2)/3.0)
  i = int(0.5*(float(x)/r - ((j+k)%2)))
  return (i,j)
