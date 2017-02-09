#!python 
# distutils: language=c++

import numpy as np
cimport numpy as np
cimport cython

from CellList import CellList
from libc.math cimport fabs  as c_fabs
from libc.math cimport ceil  as c_ceil
from libc.math cimport floor as c_floor

cdef class Box:
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
  def __init__(self,L=10,cell_grid=None):
    self.system = None

    try:
      lx = L[0]
      ly = L[1]
      lz = L[2]
    except TypeError:
      lx = ly = lz = L

    self._lx=lx
    self._ly=ly
    self._lz=lz
    self._half_lx=lx/2.0
    self._half_ly=ly/2.0
    self._half_lz=lz/2.0
    self._xlo =-lx/2.0
    self._xhi = lx/2.0
    self._ylo =-ly/2.0
    self._yhi = ly/2.0
    self._zlo =-lz/2.0
    self._zhi = lz/2.0
    if cell_grid is not None:
      self.neighbor_list = CellList(*cell_grid)
      self.neighbor_list.set_box_size(self._lx,self._ly,self._lz)
    else:
      self.neighbor_list = None

  def __str__(self):
    xyz = ( '{}:{:5.4f} '*3).format('x',self.lx,'y',self.ly,'z',self.lz)
    return '< ' + xyz + '>'

  def fit(self,positions):
    ''' Fit box to currently loaded positions in the system. '''
    maxPos = np.max(np.abs(positions),axis=0)
    self.lx=2*(maxPos[0]+1)
    self.ly=2*(maxPos[1]+1)
    self.lz=2*(maxPos[2]+1)

  def setVolume(self,vol):
    ''' Set the box lengths isotropically based on a desired volume. '''
    self.L = vol**(1.0/3.0)

  def _setLength(self,dim,length):
    if dim=='x':
      self._lx = length
      self._half_lx = length/2.0
      self._xlo = -length/2.0
      self._xhi = +length/2.0
    elif dim=='y':
      self._ly = length
      self._half_ly = length/2.0
      self._ylo = -length/2.0
      self._yhi = +length/2.0
    elif dim=='z':
      self._lz = length
      self._half_lz = length/2.0
      self._zlo = -length/2.0
      self._zhi = +length/2.0
    if self.neighbor_list is not None:
      self.neighbor_list.set_box_size(self._lx,self._ly,self._lz)

  def _setEdge(self,dim,side,value):
    if dim=='x':
      if side=='lo':
        self._xlo = value
      elif side=='hi':
        self._xhi = value
      length = (self._xhi - self._xlo)
      self._lx = length
      self._half_lx = length/2.0
    elif dim=='y':
      if side=='lo':
        self._ylo = value
      elif side=='hi':
        self._yhi = value
      length = (self._yhi - self._ylo)
      self._ly = length
      self._half_ly = length/2.0
    elif dim=='z':
      if side=='lo':
        self._zlo = value
      elif side=='hi':
        self._zhi = value
      length = (self._zhi - self._zlo)
      self._lz = length
      self._half_lz = length/2.0
    if self.neighbor_list is not None:
      self.neighbor_list.set_box_size(self._lx,self._ly,self._lz)

  ###############
  # BOX LENGTHS #
  ###############
  @property
  def L(self):
    return (self._lx,self._ly,self._lz)
  @L.setter
  def L(self,value):
    '''
    If value is a list/tuple, set the lengths independently, otherwise
    set all lengths equal to value.
    '''
    try:
      lx = value[0]
      ly = value[1]
      lz = value[2]
    except TypeError:
      lx = ly = lz = value
    self._setLength('x',lx)
    self._setLength('y',ly)
    self._setLength('z',lz)
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

  ##################
  # BOX BOUNDARIES #
  ##################
  @property
  def xlo(self):
    return self._xlo
  @xlo.setter
  def xlo(self,value):
    self._setEdge('x','lo',value)
  @property
  def ylo(self):
    return self._ylo
  @ylo.setter
  def ylo(self,value):
    self._setEdge('y','lo',value)
  @property
  def zlo(self):
    return self._zlo
  @zlo.setter
  def zlo(self,value):
    self._setEdge('z','lo',value)
  @property
  def xhi(self):
    return self._xhi
  @xhi.setter
  def xhi(self,value):
    self._setEdge('x','hi',value)
  @property
  def yhi(self):
    return self._yhi
  @yhi.setter
  def yhi(self,value):
    self._setEdge('y','hi',value)
  @property
  def zhi(self):
    return self._zhi
  @zhi.setter
  def zhi(self,value):
    self._setEdge('z','hi',value)

  @property
  def half_L(self):
    return (self._half_lx,self._half_ly,self._half_lz)
  @property
  def half_lx(self):
    return self._half_lx
  @property
  def half_ly(self):
    return self._half_ly
  @property
  def half_lz(self):
    return self._half_lz

  ###################
  # PYTHON WRAPPERS #
  ###################
  def wrap_all_positions(self,wrap_long=False):
    ''' Convenience function to wrapping all system positions back into the box. '''

    if wrap_long:
      (x,y,z),(imx,imy,imz) = self.wrap_positions_long(self.system.x,self.system.y,self.system.z)
    else:
      (x,y,z),(imx,imy,imz) = self.wrap_positions(self.system.x,self.system.y,self.system.z)
    self.system.x    = np.array(x)
    self.system.y    = np.array(y)
    self.system.z    = np.array(z)
    self.system.imx += imx
    self.system.imy += imy
    self.system.imz += imz
    self.system.reset_all_molecules()

  def wrap_positions(self,double[:] x, double[:] y, double[:] z):
    ''' Wraps *coordinate* up to **one** image distance away from the central image

    .. Warning::
      This function only wraps "once". This means that if a bead is several
      box distances away from the central box in one or more directions, this function
      would have to be called *multiple* times to wrap the position into the central
      image. Use :func:`typySim.core.cy.Box.Box.wrap_positions_long` if there is a chance 
      of this situation



    '''
    cdef long natoms = x.shape[0]
    cdef long[:] imx = np.zeros(natoms,dtype=int)
    cdef long[:] imy = np.zeros(natoms,dtype=int)
    cdef long[:] imz = np.zeros(natoms,dtype=int)
    cdef double[:] nx = x.copy()
    cdef double[:] ny = y.copy()
    cdef double[:] nz = z.copy()
    cdef long i
    for i in range(natoms):
      if nx[i] > self._xhi:
        nx[i]    -= self._lx
        imx[i]   -= 1
      elif nx[i] < self._xlo:
        nx[i]    += self._lx
        imx[i]   += 1

      if ny[i] > self._yhi:
        ny[i]    -= self._ly
        imy[i]   -= 1
      elif ny[i] < self._ylo:
        ny[i]    += self._ly
        imy[i]   += 1

      if nz[i] > self._zhi:
        nz[i]    -= self._lz
        imz[i]   -= 1
      elif nz[i] < self._zlo:
        nz[i]    += self._lz
        imz[i]   += 1
    return (nx,ny,nz),(imx,imy,imz)

  def wrap_positions_long(self,double[:] x, double[:] y, double[:] z):
    ''' Wraps *coordinate* into the central image

    .. Warning::
      This function is less efficient than :func:`self.wrap_positions` and should
      and should be avoided unless long wrapping is necessary.
    '''
    cdef long natoms = x.shape[0]
    cdef long[:] imx = np.zeros(natoms,dtype=int)
    cdef long[:] imy = np.zeros(natoms,dtype=int)
    cdef long[:] imz = np.zeros(natoms,dtype=int)
    cdef double[:] nx = x.copy()
    cdef double[:] ny = y.copy()
    cdef double[:] nz = z.copy()
    cdef long i
    cdef long image #image will be the signed image the unwrapped bead is in
    for i in range(natoms):
      if nx[i] > self._xhi:
        image     = <long> c_ceil((nx[i] - self._xhi)/self._lx)
        nx[i]    -= image*self._lx
        imx[i]   += image
      elif nx[i] < self._xlo:
        image     = <long> c_floor((nx[i] - self._xlo)/self._lx)
        nx[i]    -= image*self._lx
        imx[i]   += image

      if ny[i] > self._yhi:
        image     = <long> c_ceil((ny[i] - self._yhi)/self._ly)
        ny[i]    -= image*self._ly
        imy[i]   += image
      elif ny[i] < self._ylo:
        image     = <long> c_floor((ny[i] - self._ylo)/self._ly)
        ny[i]    -= image*self._ly
        imy[i]   += image

      if nz[i] > self._zhi:
        image     = <long> c_ceil((nz[i] - self._zhi)/self._lz)
        nz[i]    -= image*self._lz
        imz[i]   += image
      elif nz[i] < self._zlo:
        image     = <long> c_floor((nz[i] - self._zlo)/self._lz)
        nz[i]    -= image*self._lz
        imz[i]   += image

    return (nx,ny,nz),(imx,imy,imz)

  def wrap_distances(self,double[:] dx, double[:] dy, double[:] dz):
    ''' Wraps *distances* up to one box length in each direction.

    .. Warning::
      This function only wraps "once". This means that if a bead is several
      box distances away from the central box in one or more directions, this function
      would have to be called *multiple* times to wrap the position into the central
      image. 

    '''
    cdef long natoms = dx.shape[0]
    cdef double idx,idy,idz
    cdef long i
    cdef list dist = []

    for i in range(natoms):
      idx = c_fabs(dx[i])
      if idx > self._half_lx:
        idx    -= self._lx

      idy = c_fabs(dy[i])
      if idy > self._half_ly:
        idy    -= self._ly

      idz = c_fabs(dz[i])
      if idz > self._half_lz:
        idz    -= self._lz

      dist.append(((idx*idx) + (idy*idy) + (idz*idz))**(0.5))

    return dist

  def wrap_distances_long(self,double[:] dx, double[:] dy, double[:] dz):
    ''' Wraps *distances*

    .. Warning::
      This function is less efficient than :func:`self.wrap_distances` and should
      and should be avoided unless long wrapping is necessary.

    '''
    cdef long natoms = dx.shape[0]
    cdef double idx,idy,idz
    cdef long i
    cdef double image
    cdef list dist = []

    for i in range(natoms):
      idx = c_fabs(dx[i])
      if idx > self._half_lx:
        image   = c_ceil((idx - self._xhi)/self._lx)
        idx    -= image*self._lx


      idy = c_fabs(dy[i])
      if idy > self._half_ly:
        image   = c_ceil((idy - self._yhi)/self._ly)
        idy    -= image*self._ly

      idz = c_fabs(dz[i])
      if idz > self._half_lz:
        image   = c_ceil((idz - self._zhi)/self._lz)
        idz    -= image*self._lz

      dist.append(((idx*idx) + (idy*idy) + (idz*idz))**(0.5))

    return dist

  ########################
  # FAST CYTHON WRAPPERS #
  ########################
  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef double wrap_x(self,double x) nogil :
    ''' Wrap a passed *coordinate* back into the box.

    .. warning:: 
      This method **does not** update the system image flags. Use with caution. 
    '''
    if x>self._xhi:
      x -=self._lx
    elif x<self._xlo:
      x +=self._lx
    return x

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef double wrap_y(self,double y) nogil :
    ''' Wrap a passed *coordinate* back into the box.

    .. warning:: 
      This method **does not** update the system image flags. Use with caution. 
    '''
    if y>self._yhi:
      y -=self._ly
    elif y<self._ylo:
      y +=self._ly
    return y

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef double wrap_z(self,double z) nogil :
    ''' Wrap a passed *coordinate* back into the box.

    .. warning:: 
      This method **does not** update the system image flags. Use with caution. 
    '''
    if z>self._zhi:
      z -=self._lz
    elif z<self._zlo:
      z +=self._lz
    return z

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef double wrap_dx(self,double dx) nogil :
    ''' Wrap a passed *distance* back into the box.'''
    dx = c_fabs(dx)
    if dx>self._half_lx:
      dx -=self._lx
    return dx

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef double wrap_dy(self,double dy) nogil :
    ''' Wrap a passed *distance* back into the box.'''
    dy = c_fabs(dy)
    if dy>self._half_ly:
      dy -=self._ly
    return dy

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef double wrap_dz(self,double dz) nogil :
    ''' Wrap a passed *distance* back into the box.'''
    dz = c_fabs(dz)
    if dz>self._half_lz:
      dz -=self._lz
    return dz

