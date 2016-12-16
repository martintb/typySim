#!python
# distutils: language=c++

cdef class CellList(object):
  cdef bint _ready
  cdef object logger
  cdef bint central_origin
  cdef long nx,ny,nz,
  cdef public double dx,dy,dz,bx,by,bz
  cdef long ncells,ncells_1d,ncells_2d,ncells_3d
  cdef long nbeads
  cdef long[:] top
  cdef long[:] neigh
  cdef long[:] bead_cells
  cdef long[:,:] cell_neighs
  cpdef void insert_bead(self,long beadNo, double x,double y,double z) except *
  cpdef void remove_bead(self,long beadNo) except *
  cpdef void update_bead(self,long beadNo, double x,double y,double z) except *
  cdef double wrap_position(self,double x, double bx, bint central_origin) nogil
  cdef long pos2idex(self,double x, double dx, double bx, bint central_origin) nogil
  cdef long idex2cell(self,long ix,long iy,long iz) nogil
  cdef void cell2idex(self,long cell_number,long[:] ixyz) 
  cpdef void get_cell_neighbors(self,long cell_number,long[:] cell_neighs) except *
    

    
    
      







