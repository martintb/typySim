#!python 
# distutils: language=c++

cdef class Box:
  cdef object system,neighbor_list
  cdef float _lx,_ly,_lz
  cdef float _half_lx,_half_ly,_half_lz
  cdef float _xlo,_xhi
  cdef float _ylo,_yhi
  cdef float _zlo,_zhi
  cdef void wrap_position(self,double x,double y, double z)
  cdef void wrap_distance(self,double dx,double dy, double dz)
