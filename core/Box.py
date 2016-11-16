import numpy as np
from ..cy import CellList
class Box(object):
  '''
  Box utility which handles perodic wrapping and coordinate scaling.
  '''
  def __init__(self,L=-1,system=None,cell_grid=None):
    self.system = system
    self.beadVolume=0
    self._lx=L
    self._ly=L
    self._lz=L
    self._half_lx=L/2.0
    self._half_ly=L/2.0
    self._half_lz=L/2.0
    self.xlo =-L/2.0
    self.xhi = L/2.0
    self.ylo =-L/2.0
    self.yhi = L/2.0
    self.zlo =-L/2.0
    self.zhi = L/2.0
    if cell_grid is not None:
      self.cellList = CellList(*cell_grid)
    else:
      self.cellList = None
  def wrap_all_positions(self):
    self.system.positions = self.wrap_position(self.system.positions)
    self.system.reset_all_molecules()
  def wrap_position(self,vec):
    wrapped_vec = vec.copy()
    wrapped_vec = np.where(wrapped_vec>self.half_L,wrapped_vec-self.L,wrapped_vec)
    wrapped_vec = np.where(wrapped_vec<self.negative_half_L,wrapped_vec+self.L,wrapped_vec)
    return wrapped_vec
  def wrap_distance(self,dr):
    dr = np.abs(dr)
    dr = np.where(dr>self.half_L,dr-self.L,dr)
    dist = np.sqrt(np.sum(np.square(dr),axis=1))
    return dist
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
    setattr(self,'_half_l{}'.format(dim),length/2.0)
    setattr(self,'{}lo'.format(dim),-length/2.0)
    setattr(self,'{}hi'.format(dim), length/2.0)
    if self.cellList is not None:
      self.cellList.set_box_size(self._lx,self._ly,self._lz)
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
  @property
  def half_lx(self):
    return self._half_lx
  @property
  def half_ly(self):
    return self._half_ly
  @property
  def half_lz(self):
    return self._half_lz
  @property
  def half_L(self):
    return (self._half_lx,self._half_ly,self._half_lz)
  @property
  def negative_half_L(self):
    return (-self._half_lx,-self._half_ly,-self._half_lz)

