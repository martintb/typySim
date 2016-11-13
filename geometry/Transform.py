import numpy as np
from linalg import rotation_matrix

class Transform(object):
  def __init__(self):
    self.box = None
    self.anchor = None #center point for rotations/translations
  def translate(self,new_anchor_pos,positions):
    if self.anchor is None:
      raise ValueError('Cannot apply transform without an anchor!')
    v = new_anchor_pos-self.anchor
    new_positions = positions + v
    return new_positions
  def rotate(self,v1,v2,positions):
    R = rotation_matrix(v1,v2)
    # this assumes positions is an Nx3 array
    new_positions = np.dot(positions[np.newaxis,:,:],R)
    return new_positions




