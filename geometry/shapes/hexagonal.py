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
  for k in range(nz):                      
    for j in range(ny):                    
      for i in range(nx):                  
        x = (2*i + ((j+k)%2))*radius
        y = sqrt(3)*(j+(k%2)/3.0)*radius
        z = 2.0/3.0 * sqrt(6) * k * radius
        positions.append([x,y,z])
        if k==0:
          types.append(0)
        elif k==(nz-1):
          types.append(1)
        else:
          types.append(2)
  positions = np.array(positions)
  if shift:
    positions -= np.average(positions,axis=0)
  molData ={}
  molData['positions'] = positions
  molData['types'] = types
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
