
cdef inline double norm(double[:])
cdef inline double normalize(double[:])
cdef inline double dot(double[:],double[:])
cdef inline void cross(double[:],double[:],double[:])
cdef inline double cross_magnitude(double[:],double[:])
cdef inline void subtract(double[:],double[:],double[:],double[:])
cdef inline double periodic_boundary_wrap(double, double)
