from math import sqrt
import numpy as np

def line(num,diameter,bonds=None,angles=None):
  positions = [[0.,0.,i*diameter] for i in range(num)]
  if bonds is not None:
    bonds = [[i,j] for i,j in zip(range(num-1),range(1,num))]
  if angles is not None:
    angles = [[i,j,k] for i,j,k in zip(range(num-2),range(1,num-1),range(2,num))]
  return positions

def hexagonal_surface(nx,ny,nz,diameter,gen_types=False,shift=True):
  radius = diameter/2.0                    
  positions = []                           
  if gen_types is True:
    types = []
  for k in range(nz):                      
    for j in range(ny):                    
      for i in range(nx):                  
        x = (2*i + ((j+k)%2))*radius
        y = sqrt(3)*(j+(k%2)/3.0)*radius
        z = 2.0/3.0 * sqrt(6) * k * radius
        positions.append([x,y,z])
        if gen_types is True:
          if k==0:
            types.append(0)
          elif k==(nz-1):
            types.append(1)
          else:
            types.append(2)
  positions = np.array(positions)
  if shift:
    positions -= np.average(positions,axis=0)
  result ={}
  result['positions'] = positions
  if gen_types is True:
    result['types'] = types
  return result
