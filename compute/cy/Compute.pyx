#!python
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False


cdef class Compute:
  def __init__(self):
    self.system = None
    self._name = "BaseCompute"
    self.frame_skip = -1
    self.block_size = -1
    self.num_blocks = -1
    self.values = []
  @property
  def name(self):
    return self._name
  def set_calc_times(self,frame_skip,block_size,num_frames):
    self.frame_skip = frame_skip
    self.block_size = block_size
    self.num_frames = num_frames
  def compute(self):
    return NotImplementedError('This compute has not defined `compute`!')

