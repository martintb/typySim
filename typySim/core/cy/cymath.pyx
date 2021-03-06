#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
from libc.math cimport sqrt as c_sqrt

from typySim.core.cy.Box cimport * 

cdef inline double norm(double[:] vec):
  return c_sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])

cdef inline double normalize(double[:] vec):
  cdef double normVal = c_sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]) 
  vec[0] /= normVal
  vec[1] /= normVal
  vec[2] /= normVal

cdef inline double dot(double[:] u, double[:] v): 
  return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

cdef inline void cross(double[:] u, double[:] v, double[:] uxv):
  uxv[0] = u[1]*v[2] - u[2]*v[1]
  uxv[1] = u[2]*v[0] - u[0]*v[2]
  uxv[2] = u[0]*v[1] - u[1]*v[0]

cdef inline double cross_magnitude(double[:] u, double[:] v):
  return c_sqrt(\
                (u[1]*v[2] - u[2]*v[1]) * (u[1]*v[2] - u[2]*v[1]) +   \
                (u[2]*v[0] - u[0]*v[2]) * (u[2]*v[0] - u[0]*v[2]) +   \
                (u[0]*v[1] - u[1]*v[0]) * (u[0]*v[1] - u[1]*v[0])     \
               )

cdef inline void subtract(double[:] u, double[:] v, double[:] uv, double[:] BOXL):
  uv[0] = periodic_boundary_wrap(u[0]-v[0], BOXL[0])
  uv[1] = periodic_boundary_wrap(u[1]-v[1], BOXL[1])
  uv[2] = periodic_boundary_wrap(u[2]-v[2], BOXL[2])

cdef inline double periodic_boundary_wrap(double dr, double BOXL):
  cdef double B2 =  BOXL/2.0
  if dr>B2:
    dr=dr-BOXL
  elif dr<-B2:
    dr=dr+BOXL
  return dr

# cdef long binary_search(long value, long[:] array) nogil:
#   ''' 
#   Returns index of value in array if value is found
#   Return -1 is value is not found
#   Warning: Array must be pre-sorted for this to work! 
#   '''
#   cdef long lo = 0
#   cdef long hi = array.shape[0]-1
#   cdef long mid
#   while lo<=hi:
#     mid = lo + (hi-lo)/2
#     if array[mid] == value:
#       return mid
#     elif array[mid] < value:
#       lo = mid + 1
#     else:
#       hi = mid - 1
#   #value not found!
#   return -1
# 
# cpdef n_closest(long n, double x1,double y1,double z1, double[:] x2,double[:] y2,double[:] z2, Box box):
#   ''' 
#   Returns 
#   '''
#   cdef long[:] index_out = np.full(n,-1,dtype=np.float)
#   cdef double[:] dist_out = np.full(n,-1,dtype=np.int)
# 
#   cdef long natoms = x2.shape[0]
#   cdef long i,j
# 
#   return np.array(index_out),np.array(dist_out)
