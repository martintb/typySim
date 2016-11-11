import numpy as np
class Box(object):
  def __init__(self,L=-1):
    self.beadVolume=0
    self._lx=L
    self._ly=L
    self._lz=L
    self.xlo =-L/2.0
    self.xhi = L/2.0
    self.ylo =-L/2.0
    self.yhi = L/2.0
    self.zlo =-L/2.0
    self.zhi = L/2.0
  def fit(self,positions,iso=True):
    if iso:
      maxPos = np.max(np.abs(positions))
      self.L=2*(maxPos+1)
    else:
      maxPos = np.max(np.abs(positions),axis=0)
      self.lx=2*(maxPos[0]+1)
      self.ly=2*(maxPos[1]+1)
      self.lz=2*(maxPos[2]+1)
  def setVolume(self,vol):
    self.L = vol**(1.0/3.0)
  def setVolumeFraction(self,vFrac,volumeOffset=None):
    vol = self.beadVolume
    if volumeOffset is not None:
      vol += volumeOffset
    self.L = (vol/vFrac)**(1.0/3.0)
  def setLength(self,dim,length):
    setattr(self,'_l{}'.format(dim),length)
    setattr(self,'{}lo'.format(dim),-length/2.0)
    setattr(self,'{}hi'.format(dim), length/2.0)
  def __str__(self):
    xyz = ( '{}:{:5.4f} '*3).format('x',self.lx,'y',self.ly,'z',self.lz)
    return '< ' + xyz + '>'
  @property
  def L(self):
    return (self._lx,self._ly,self._lz)
  @L.setter
  def L(self,value):
    self.setLength('x',value)
    self.setLength('y',value)
    self.setLength('z',value)
  @property
  def lx(self):
    return self._lx
  @lx.setter
  def lx(self,value):
    self.setLength('x',value)
  @property
  def ly(self):
    return self._ly
  @ly.setter
  def ly(self,value):
    self.setLength('y',value)
  @property
  def lz(self):
    return self._lz
  @lz.setter
  def lz(self,value):
    self.setLength('z',value)
