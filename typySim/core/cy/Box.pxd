#!python 
# distutils: language=c++

cdef class Box(object):
  cdef public object system
  cdef public object neighbor_list
  cdef float _lx,_ly,_lz
  cdef float _half_lx,_half_ly,_half_lz
  cdef float _xlo,_xhi
  cdef float _ylo,_yhi
  cdef float _zlo,_zhi
  cdef double wrap_x(self,double) nogil
  cdef double wrap_y(self,double) nogil
  cdef double wrap_z(self,double) nogil
  cdef double wrap_dx(self,double) nogil
  cdef double wrap_dy(self,double) nogil
  cdef double wrap_dz(self,double) nogil
