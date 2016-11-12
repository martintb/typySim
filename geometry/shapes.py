from math import sqrt
import numpy as np

def line(length,diameter):
  molData = {}
  molData['positions'] = [[0.,0.,i*diameter] for i in range(length)]
  molData['bonds'] = [[i,j] for i,j in zip(range(length-1),range(1,length))]
  molData['angles'] = [[i,j,k] for i,j,k in zip(range(length-2),range(1,length-1),range(2,length))]
  return molData

def hexagonal_surface(nx,ny,nz,diameter,shift=True):
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
