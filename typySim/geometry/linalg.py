import numpy as np


def normalize(v,axis=None):
  v /= np.linalg.norm(v,axis=axis)
  return v

def rescale(v,scale,axis=None):
  v = np.linalg.norm(v,axis=axis)*scale
  return v

class Quaternion(object):
    def __init__(self,qvector=None,axis=None,angle=None,vec_from=None,vec_to=None):
        if qvector is not None:
            self.qvector = qvector
        elif (axis is not None) and (angle is not None):
            self.set_rotation(axis=axis,angle=angle)
        elif (vec_from is not None) and (vec_to is not None):
            axis,angle = self.get_axis_angle(vec_from=vec_from,vec_to=vec_to)
            self.set_rotation(axis=axis,angle=angle)
        else:
            self.qvector = (0,0,0,0)
    def get_axis_angle(self,vec_from,vec_to):
        vec_from = vec_from/np.linalg.norm(vec_from)
        vec_to   = vec_to/np.linalg.norm(vec_to)

        axis     = np.cross(vec_from,vec_to)
        angle    = np.arccos(np.dot(vec_from,vec_to))
        return axis,angle
    def set_rotation(self,axis,angle):
        u = axis/np.linalg.norm(axis)
        theta = angle/2.0
        w = np.cos(theta)
        x = u[0]*np.sin(theta)
        y = u[1]*np.sin(theta)
        z = u[2]*np.sin(theta)
        self.qvector = (w,x,y,z)
    @property
    def inverse_qvector(self):
        return (self.iw,self.ix,self.iy,self.iz)
    @property
    def qvector(self):
        return (self.w,self.x,self.y,self.z)
    @qvector.setter
    def qvector(self,qvector):
        self.w,self.x,self.y,self.z = qvector
        self.iw = +1.0*self.w
        self.ix = -1.0*self.x
        self.iy = -1.0*self.y
        self.iz = -1.0*self.z
    def __mul__(self,other):
        q1 = self.qvector
        q2 = other.qvector
        q3 = self.hamilton_product(q1,q2)
        return Quaternion(qvector=q3)
    def hamilton_product(self,q1,q2):
        w1,x1,y1,z1 = q1
        w2,x2,y2,z2 = q2
        w = w1*w2 - x1*x2 - y1*y2 - z1*z2
        x = w1*x2 + x1*w2 + y1*z2 - z1*y2
        y = w1*y2 - x1*z2 + y1*w2 + z1*x2
        z = w1*z2 + x1*y2 - y1*x2 + z1*w2
        return w,x,y,z
    def rotate(self,vector,rotation_axis=None,rotation_angle=None):
        if (rotation_axis is not None) and (rotation_angle is not None):
            self.set_rotation(rotation_axis,rotation_angle)
        q1 = self.qvector
        q2 = (0,) + tuple(vector)
        q3 = self.inverse_qvector
        q4 = self.hamilton_product(self.hamilton_product(q1,q2),q3)
        return np.array(q4[1:])
        
