from hexagonal import position2Index,index2Position
import numpy as np
# import ipdb; ist = ipdb.set_trace

def overlap(r1,r2,contact_dist=1.0):
  if r1.size!=3:
    raise ValueError('This function only works if r1 is a single coordinate')

  allDist =  np.sqrt(np.sum(np.square(r1 - r2),axis=1))

  if np.any(allDist<(contact_dist-0.0000001)):
    return True
  else:
    return False

def overlapxy(r1xyz,r2xyz,contact_dist=1.0):
  if r1xyz.size!=3:
    raise ValueError('This function only works if r1 is a single coordinate')

  r1 = r1xyz[:,:2]
  r2 = r2xyz[:,:2]

  allDist =  np.sqrt(np.sum(np.square(r1 - r2),axis=1))

  if np.any(allDist<(contact_dist-0.0000001)):
    return True
  else:
    return False

def positions(natoms, d=1.0):
  a=d; b=0.5*d; c=b*(3.**(0.5))
  baseHex = [[a,0.,0.],[b,c,0.],[-b,c,0.],[-a,0.,0.],[-b,-c,0.],[b,-c,0.]]

  stepDir=0
  pos=np.array([[0.,0.,0.]])
  for i in range(natoms-1):

    r1=np.array([pos[i]])
    hexDir=(stepDir+2)%6
    r2=np.array([pos[i]+baseHex[hexDir]])
    hexDir=(stepDir+1)%6
    r3=np.array([pos[i]+baseHex[hexDir]])
    hexDir=(stepDir)%6
    r4=np.array([pos[i]+baseHex[hexDir]])

    if not overlapxy(r2,pos):
      nextPos = r2
      stepDir = (stepDir+2)%6
    elif not overlapxy(r3,pos):
      nextPos = r3
      stepDir=(stepDir+1)%6
    else:
      nextPos = r4
    i,j,_ = position2Index(nextPos[0][0],nextPos[0][1],r=d/2.0,alternate_z=True)
    nextPos[0][2] = index2Position(abs(i),abs(j),r=d/2.0,alternate_z=True)[2]
    pos = np.append(pos,nextPos,axis=0)
  return pos

