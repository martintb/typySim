import os

class Writer(object):
  '''
  Base class for all writers
  '''
  def __init__(self,fileBase):
    self.name = 'BaseWriter'
    self.file = None
    self.fileBase = fileBase
  def write_frame(*args,**kwargs):
    raise NotImplementedError('write_frame shouldn\'t be called on the parent class!')
  
