#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

cimport numpy as np
import numpy as np

from typySim.core.cy.Box cimport * 

from libc.math cimport sqrt as c_sqrt

cdef long binary_search(long value, long[:] array) nogil:
  ''' 
  Returns index of value in array if value is found
  Return -1 is value is not found
  Warning: Array must be pre-sorted for this to work! 
  '''
  cdef long lo = 0
  cdef long hi = array.shape[0]-1
  cdef long mid
  while lo<=hi:
    mid = lo + (hi-lo)/2
    if array[mid] == value:
      return mid
    elif array[mid] < value:
      lo = mid + 1
    else:
      hi = mid - 1
  #value not found!
  return -1

cpdef n_closest(long n, double x1,double y1,double z1, double[:] x2,double[:] y2,double[:] z2, Box box):
  ''' 
  Returns 
  '''
  cdef long[:] index_out = np.full(n,-1,dtype=np.int)
  cdef double[:] dist_out = np.full(n,1e6,dtype=np.float)

  cdef long natoms = x2.shape[0]
  cdef long i,j
  cdef double dx,dy,dz,dist

  for i in range(natoms):
    dx = x1 - x2[i]
    dy = y1 - y2[i]
    dz = z1 - z2[i]

    dx = box.wrap_dx(dx)
    dy = box.wrap_dy(dy)
    dz = box.wrap_dz(dz)

    dist = c_sqrt(dx*dx + dy*dy + dz*dz)
    for j in range(n):
      if dist<dist_out[j]:
        dist_out[j] = dist
        index_out[j] = i
        break

  return np.array(index_out),np.array(dist_out)
