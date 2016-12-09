#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport fabs as c_fabs
from libc.math cimport sqrt as c_sqrt
import logging

intType = np.long
floatType = np.float32
doubleType = np.float
ctypedef np.int_t cIntType
ctypedef np.float32_t cFloatType
ctypedef np.float_t cDoubleType


cdef class CellList:
  '''
  Implementation of the doubly-linked list based cell-list algorithm


  cell indicies example
  ---------------------------------
  |(0,0,0)|(0,1,0)|(0,2,0)|(0,3,0)|
  ---------------------------------
  |(1,0,0)|(1,1,0)|(1,2,0)|(1,3,0)|
  ---------------------------------
  |(2,0,0)|(2,1,0)|(2,2,0)|(2,3,0)|
  ---------------------------------
  |(3,0,0)|(3,1,0)|(3,2,0)|(3,3,0)|
  ---------------------------------
  |(4,0,0)|(4,1,0)|(4,2,0)|(4,3,0)|
  ---------------------------------

  cell number example
  ---------------------------------
  |  (0)  |  (1)  |  (2)  |  (3)  |
  ---------------------------------
  |  (4)  |  (5)  |  (6)  |  (7)  |
  ---------------------------------
  |  (8)  |  (9)  | (10)  | (11)  |
  ---------------------------------
  | (12)  | (13)  | (14)  | (15)  |
  ---------------------------------
  | (16)  | (17)  | (18)  | (19)  |
  ---------------------------------

  Terminology
  -----------
    cell number: single integer which identifies a cell
    cell index:  triple of integers which identifies a cell
  '''
  def __init__(self,nx,ny,nz):
    if (nx<3) or (ny<3) or (nz<3):
      raise ValueError('Need at least 3 divisions in each direction for CellList to work!')
    self._ready=False
    self.bx = -1
    self.by = -1
    self.bz = -1
    self.dx = -1
    self.dy = -1
    self.dz = -1
    self.nbeads = 0
    self.nx = nx
    self.ny = ny
    self.nz = nz
    self.ncells = nx*ny*nz
    self.ncells_1d = nx
    self.ncells_2d = nx*ny
    self.ncells_3d = nx*ny*nz
    self.logger = logging.getLogger(__name__)
    self.logger.debug('==> Cython CellList Object Created:')
    self.logger.debug('--> Cell Grid Length: {:d} {:d} {:d}'.format(self.nx,self.ny,self.nz))
    self.logger.debug('--> Number of Cells: {:d}'.format(self.ncells))
    self.logger.debug('--> Number of Cells 1D,2D,3D: {:d} {:d} {:d}'.format(self.ncells_1d,self.ncells_2d,self.ncells_3d))
    self.create_cell_neighbors()
  @property
  def ready(self):
    return self._ready
  def create_cell_neighbors(self):
    '''
    Create a list of each cells immediate "neighbors"

    For performance reasons, we pre-allocate a the list of neighbors
    of each cell, so that we can quickly iterate over it when calculating
    neighbor distances. 

    '''
    self.logger.debug('.:: Initializing cell neighbors...')
    cdef long cellno,m
    cdef long sx,sy,sz #shift xyz
    cdef long ix,iy,iz #index xyz
    cdef long nix,niy,niz #neigh index xyz
    cdef long[:] ixyz = np.full(3,-1,dtype=intType)
    self.cell_neighs = np.full((self.ncells,27),-1,dtype=intType)
    for cellno in range(self.ncells):
      self.cell2idex(cellno,ixyz)
      ix = ixyz[0]
      iy = ixyz[1]
      iz = ixyz[2]
      m=0
      for sx in [-1,0,1]:
        for sy in [-1,0,1]:
          for sz in [-1,0,1]:
            nix = ix + sx
            if nix>=self.nx:
              nix  = 0
            elif nix<0:
              nix = self.nx-1
            niy = iy + sy
            if niy>=self.ny:
              niy  = 0
            elif niy<0:
              niy = self.ny-1
            niz = iz + sz
            if niz>=self.nz:
              niz  = 0
            elif niz<0:
              niz = self.nz-1
            self.cell_neighs[cellno,m] = self.idex2cell(nix,niy,niz)
            # print cellno,m,ix,iy,iz,nix,niy,niz,self.cell_neighs[cellno,m]
            m+=1
    self.logger.debug('==> Done initializing cell neighbors')
  def reset_nlist(self,nbeads):
    '''
    Reset/resize the neighborlist, overwriting any current contents.
    '''
    self.logger.debug('==> Setting up CellList..')
    self.logger.debug('--> Number of beads: {:d}'.format(nbeads))
    self.logger.debug('--> Box Size: {:3.2f} {:3.2f} {:3.2f}'.format(self.bx,self.by,self.bz))
    self.logger.debug('--> Cell Grid Spacing: {:3.2f} {:3.2f} {:3.2f}'.format(self.dx,self.dy,self.dz))

    self.nbeads = nbeads
    self.top = np.full(self.ncells,-1, dtype=intType)
    self.neigh = np.full(self.nbeads,-1, dtype=intType)
    self.bead_cells = np.full(self.nbeads,-1, dtype=intType)
  def expand(self,num):
    '''
    Expand the neighborlist by :param num:, keeping the current contents in place
    '''
    self.logger.debug('==> Resizing CellList by {}'.format(num))
    self.neigh = np.append(self.neigh,[-1]*num)
    self.bead_cells = np.append(self.bead_cells,[-1]*num)
    self.nbeads+=num
  def shrink(self,indices):
    '''
    Shrink the neighborlist by removing :param indices:,
    '''
    self.logger.debug('==> Resizing CellList by removing {}'.format(indices))
    self.neigh = np.delete(self.neigh,indices)
    self.bead_cells = np.delete(self.bead_cells,indices)
    self.nbeads -= len(indices)
  def get_nlist(self):
    '''
    Pass back all neighbor_list information to user as a dictionary.
    '''
    outdict = {}
    outdict['top'] = np.array(self.top)
    outdict['neigh'] = np.array(self.neigh)
    outdict['cell_neighs'] = np.array(self.cell_neighs)
    outdict['bead_cells'] = np.array(self.bead_cells)
    outdict['nx'] = self.nx
    outdict['ny'] = self.ny
    outdict['nz'] = self.nz
    outdict['dx'] = self.dx
    outdict['dy'] = self.dy
    outdict['dz'] = self.dz
    outdict['bx'] = self.bx
    outdict['by'] = self.by
    outdict['bz'] = self.bz
    outdict['ncells'] = self.ncells
    outdict['nbeads'] = self.nbeads
    return outdict
  def set_box_size(self,double bx,double by,double bz):
    '''
    Set/Resize box and cell sizes
    '''
    self.bx = bx
    self.by = by
    self.bz = bz
    self.dx = self.bx/self.nx
    self.dy = self.by/self.ny
    self.dz = self.bz/self.nz
  def build_nlist(self,double[:] x, double[:] y, double[:] z,bint central_origin):
    '''
    Build a new neighbor_list from scratch

    To call this function, you must have already called set_box_size
    '''
    if self.bx == self.by == self.bz == -1:
      raise ValueError('This CellList doesn\'t have a box size set for it!')
    self.logger.debug('==> Creating new nlist from scratch')
    self.central_origin=central_origin
    # self.set_box_size(box[0],box[1],box[2])
    self.reset_nlist(x.shape[0])
    cdef double xx,yy,zz
    for beadNo in range(self.nbeads):
      xx = x[beadNo]
      yy = y[beadNo]
      zz = z[beadNo]
      self.insert_bead(beadNo,xx,yy,zz)
    self._ready=True
  cpdef void insert_bead(self,long beadNo, double x,double y,double z) except *:
    '''
    Insert bead into neighbor list. If neccessary, resizes the nlists.

    If beadNo<0, the function assumes you want to append a new position
    to the array. 
    '''

    cdef long ix,iy,iz,cellNo,
    cdef long oldTopBead,newTopBeadNo
    cdef long thisNeighNo,nextNeighNo

    if beadNo<0:
      beadNo = self.nbeads

    if beadNo==self.nbeads:
      self.expand(1)
    elif beadNo>self.nbeads:
      raise ValueError('To append new bead to list, beadNo must be equal to nbeads')

    ix = self.pos2idex(x,self.dx,self.bx)
    iy = self.pos2idex(y,self.dy,self.by)
    iz = self.pos2idex(z,self.dz,self.bz)
    cellNo = self.idex2cell(ix,iy,iz)
    self.bead_cells[beadNo] = cellNo

    if beadNo>self.top[cellNo]:
      # This is the new "top" bead!
      newTopBeadNo = beadNo
      oldTopBeadNo = self.top[cellNo]
      self.top[cellNo] = newTopBeadNo
      self.neigh[beadNo] = oldTopBeadNo
    else:
      # This bead belongs in the middle of the list
      # We will iterate until we find the correct location
      thisNeighNo = self.top[cellNo]
      while thisNeighNo!=-1:
        nextNeighNo = self.neigh[thisNeighNo]
        if beadNo>nextNeighNo:
          self.neigh[thisNeighNo] = beadNo
          self.neigh[beadNo] = nextNeighNo
          break
          #thisNeighNo=-1 #break loop
        else:
          thisNeighNo = nextNeighNo
  cpdef void remove_bead(self,long beadNo) except *:
    '''
    Removes a bead's location in the linked neighbor list
    '''

    cdef long ix,iy,iz,cellNo,
    cdef long oldTopBead,newTopBeadNo
    cdef long thisNeighNo,nextNeighNo

    if beadNo<0:
      beadNo = (self.nbeads-1)

    ###############################
    # remove old location in list #
    ###############################
    cellNo = self.bead_cells[beadNo]
    if beadNo==self.top[cellNo]:
      self.top[cellNo] = self.neigh[beadNo]
    else:
      # This bead is in the middle of the list. Must iterate
      thisNeighNo = self.top[cellNo]
      while thisNeighNo!=-1:
        nextNeighNo = self.neigh[thisNeighNo]
        if beadNo==nextNeighNo:
          self.neigh[thisNeighNo] = self.neigh[beadNo] #bypass this bead's node
          break
        else:
          thisNeighNo = nextNeighNo
    self.neigh[beadNo] = -1 
  cpdef void update_bead(self,long beadNo, double x,double y,double z) except *:
    '''
    Updates a bead's location in the cell list. 
    '''
    # remove old location in list 
    self.remove_bead(beadNo)

    # add new location to list 
    self.insert_bead(beadNo,x,y,z)
  cdef long pos2idex(self,double x, double dx, double bx) nogil:
    '''
    Convert coordinate to cell index
    '''
    cdef double shift
    cdef long idex
    if self.central_origin:
      shift = bx/2.0
    else:
      shift =0
    idex = <long> ((x+shift)/(dx))
    return idex
  cdef long idex2cell(self,long ix,long iy,long iz) nogil:
    '''
    Convert cell index triplet to cell number
    '''
    cdef long cellno
    cdef long n1d = self.ncells_1d
    cdef long n2d = self.ncells_2d
    cellno = <long>(ix + iy*n1d + iz*n2d)
    return cellno
  cdef void cell2idex(self,long cell_number,long[:] ixyz):
    '''
    Convert cell number to cell index triplet
    '''
    cdef long n1d = self.ncells_1d
    cdef long n2d = self.ncells_2d
    ixyz[2] = (long) (cell_number/n2d)
    ixyz[1] = (long) (cell_number- ixyz[2]*n2d)/n1d
    ixyz[0] = (long) (cell_number - ixyz[1]*n1d -ixyz[2]*n2d)
  cpdef void get_cell_neighbors(self,long cell_number,long[:] cell_neighs) except *:
    '''
    Get the cell numbers of the neighbors of a particular cell
    '''
    if cell_neighs.shape[0]!=27:
      raise TypeError('get_cell_neighbors:argument 2 must have length 27 and be of a ndarray of type long!')
    if cell_number>=self.ncells:
      raise TypeError('get_cell_neighbors:requested cell_number is higher than number of cells!')
    cdef long i
    for i in range(27):
      cell_neighs[i] = self.cell_neighs[cell_number,i]
  def get_neighbors_by_cell(self,long cellNo):
    '''
    Returns all indexes that are neighbors of requested cell (including in requested cell)
    '''
    cdef long currCell
    cdef long thisNeighNo
    cdef long i
    cdef list neighs = []
    for i in range(27):
      currCell = self.cell_neighs[cellNo,i]
      thisNeighNo = self.top[currCell]
      while thisNeighNo!=-1:
        neighs.append(thisNeighNo)
        thisNeighNo = self.neigh[thisNeighNo]
    return neighs
  def get_neighbors_by_pos(self,double x, double y, double z):
    '''
    Given a position, returns all bead indices of all neighbors
    '''
    cdef long ix,iy,iz,cellNo
    ix = self.pos2idex(x,self.dx,self.bx)
    iy = self.pos2idex(y,self.dy,self.by)
    iz = self.pos2idex(z,self.dz,self.bz)
    cellNo = self.idex2cell(ix,iy,iz)
    neighs = self.get_neighbors_by_cell(cellNo)
    return neighs
  def get_neighbors_by_index(self,long bead_index):
    '''
    Given a bead index, returns all bead indices of all neighbors
    '''
    cdef long cellNo = self.bead_cells[bead_index]
    neighs = self.get_neighbors_by_cell(cellNo)
    return neighs
  def calc_neighbor_dists(self,long bead_index,double[:] x,double[:] y,double[:] z):
    '''
    Given a bead index, returns all bead indices of all neighbors
    '''
    cdef long cellNo = self.bead_cells[bead_index]
    cdef long currCell
    cdef long thisNeighNo
    cdef long i
    cdef double x1,y1,z1
    cdef double dx,dy,dz
    cdef double Bx2 = self.bx/2.0
    cdef double By2 = self.by/2.0
    cdef double Bz2 = self.bz/2.0
    cdef double Bx = self.bx
    cdef double By = self.by
    cdef double Bz = self.bz
    cdef double rsq,r
    x1 = x[bead_index]
    y1 = y[bead_index]
    z1 = z[bead_index]
    dists = []
    indexes = []
    for i in range(27):
      currCell = self.cell_neighs[cellNo,i]
      thisNeighNo = self.top[currCell]
      while thisNeighNo!=-1:
        if thisNeighNo!=bead_index:
          dx = c_fabs (x1 - x[thisNeighNo])
          dy = c_fabs (y1 - y[thisNeighNo])
          dz = c_fabs (z1 - z[thisNeighNo])
    
          if dx>Bx2:
            dx=dx-Bx
          if dy>By2:
            dy=dy-By
          if dz>Bz2:
            dz=dz-Bz
    
          r = c_sqrt(dx*dx + dy*dy + dz*dz)
          dists.append(r)
          indexes.append(thisNeighNo)
        thisNeighNo = self.neigh[thisNeighNo]
    return np.array(dists),np.array(indexes)
    

    
    
      







