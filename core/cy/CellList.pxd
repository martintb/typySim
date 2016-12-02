#!python
# distutils: language=c++

cdef class CellList(object):
  cdef object logger
  cdef bint central_origin
  cdef long nx,ny,nz,
  cdef double dx,dy,dz,bx,by,bz
  cdef long ncells,ncells_1d,ncells_2d,ncells_3d
  cdef long nbeads
  cdef long[:] top
  cdef long[:] neigh
  cdef long[:] bead_cells
  cdef long[:,:] cell_neighs
  cpdef void insert_bead(self,long beadNo, double x,double y,double z)
  cpdef void remove_bead(self,long beadNo)
  cpdef void update_bead(self,long beadNo, double x,double y,double z)
  cpdef long pos2idex(self,double x, double dx, double bx)
  cpdef long idex2cell(self,long ix,long iy,long iz)
  cpdef void cell2idex(self,long cell_number,long[:] ixyz)
  cpdef void get_cell_neighbors(self,long cell_number,long[:] cell_neighs)
    

    
    
      







