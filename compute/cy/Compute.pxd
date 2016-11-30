#!python
# distutils: language=c++

cdef class Compute:
  cdef object system
  cdef unicode name
  cdef long frame_skip
  cdef long block_size
  cdef long num_frames
