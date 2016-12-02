#!python
# distutils: language=c++

cdef class Compute:
  cdef object system
  cdef str _name
  cdef long frame_skip
  cdef long block_size
  cdef long num_blocks
