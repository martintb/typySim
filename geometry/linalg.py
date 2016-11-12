import numpy as np


def normalize(v,axis=None):
  v /= np.linalg.norm(v,axis=axis)

def rotation_matrix(v1,v2):
  R = np.empty((3,3))
  R[:,0] = normalize(v1)
  R[:,2] = normalize(np.cross(v1,v2))
  R[:,1] = normalize(np.cross(R[:,2],v1))
  return np.array(R)
  
