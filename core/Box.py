import numpy as np
from typySim.core.cy import CellList

class Box(object):
  ''' 
  General box service which handles perodic wrapping and coordinate scaling.

  The goal of this service is to make all box length related operations 
  including resizing/rescaling, coordinate wrapping/unwrapping, and distance
  calculations sane and error-safe for the user. 

  Attributes
  ----------
  lx,ly,lz,L : float
      These are python parameters will explicitly defined setters. 
      ( e.g. setting L will automatically set lx, ly, and lz.)
  xlo,xhi,ylo,yhi,zlo,zhi : float
      Box-edge locations which are automatically calculated via the box edge
      setters. Currently, this class assumes a central-origin. Box will likely
      be extended to handle non-central origin systems in the future.

  Parameters
  ----------
  system : object
      Reference to parent system. Used for auto-wrapping and position
  cell_grid : int tuple of size, optional
      Definition of the number of divisions to use in the x, y, and z directions
      when creating a Cell List based neighorborlist object. After creation, the 
      Cell List is updated automatically via the box-edge setters. If cell_grid
      is not specified, no neighborlist is created. 

  '''
  def __init__(self,L=-1,cell_grid=None):
    self.system = None
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
      self.neighbor_list = CellList(*cell_grid)
    else:
      self.neighbor_list = None
  def wrap_all_positions(self):
    ''' Convenience function to wrapping all system positions back into the box. '''
    self.system.positions = self.wrap_position(self.system.positions)
    self.system.reset_all_molecules()
  def wrap_position(self,vec):
    ''' Wrap all passed positions back into the box.'''
    wrapped_vec = vec.copy()
    wrapped_vec = np.where(wrapped_vec>self.half_L,wrapped_vec-self.L,wrapped_vec)
    wrapped_vec = np.where(wrapped_vec<self.negative_half_L,wrapped_vec+self.L,wrapped_vec)
    return wrapped_vec
  def wrap_distance(self,dr):
    ''' Wrap all passed distances via the minimum image convention. '''
    dr = np.abs(dr)
    dr = np.where(dr>self.half_L,dr-self.L,dr)
    dist = np.sqrt(np.sum(np.square(dr),axis=1))
    return dist
  def fit(self):
    ''' Fit box to currently loaded positions in the system. '''
    maxPos = np.max(np.abs(positions),axis=0)
    self.lx=2*(maxPos[0]+1)
    self.ly=2*(maxPos[1]+1)
    self.lz=2*(maxPos[2]+1)
  def setVolume(self,vol):
    ''' Set the box lengths isotropically based on a desired volume. '''
    self.L = vol**(1.0/3.0)
  def _setLength(self,dim,length):
    setattr(self,'_l{}'.format(dim),length)
    setattr(self,'_half_l{}'.format(dim),length/2.0)
    setattr(self,'{}lo'.format(dim),-length/2.0)
    setattr(self,'{}hi'.format(dim), length/2.0)
    if self.neighbor_list is not None:
      self.neighbor_list.set_box_size(self._lx,self._ly,self._lz)
  def __str__(self):
    xyz = ( '{}:{:5.4f} '*3).format('x',self.lx,'y',self.ly,'z',self.lz)
    return '< ' + xyz + '>'
  @property
  def L(self):
    return (self._lx,self._ly,self._lz)
  @L.setter
  def L(self,value):
    self._setLength('x',value)
    self._setLength('y',value)
    self._setLength('z',value)
  @property
  def lx(self):
    return self._lx
  @lx.setter
  def lx(self,value):
    self._setLength('x',value)
  @property
  def ly(self):
    return self._ly
  @ly.setter
  def ly(self,value):
    self._setLength('y',value)
  @property
  def lz(self):
    return self._lz
  @lz.setter
  def lz(self,value):
    self._setLength('z',value)
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

