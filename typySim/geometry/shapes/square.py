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
        x = 2*i*radius
        y = 2*j*radius
        z = 2*k*radius
        positions.append([x,y,z])
        if k==0:
          types.append(topType)
        elif k==(nz-1):
          types.append(bottomType)
        else:
          types.append(middleType)
  positions = np.array(positions)
  if shift:
    positions -= np.average(positions,axis=0)
  molData ={}
  molData['x'] = positions[:,0]
  molData['y'] = positions[:,1]
  molData['z'] = positions[:,2]
  molData['types'] = types
  return molData
