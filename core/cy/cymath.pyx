#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
from libc.math cimport sqrt as c_sqrt

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
