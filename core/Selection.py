import numpy as np

class Selection(object):
  def __init__(self,select_from=None,select_vals=None):
    if (select_from is not None) and (select_vals is not None):
      self.select(select_from,select_vals)
    else:
      self.mask = None
      self.indices = None
      self.select_from = None
      self.select_vals = None
  def select(self,select_from,select_vals):
    self.select_from = select_from
    self.select_vals = select_vals
    self.mask = np.in1d(select_from,select_vals)
    self.indices = np.where(self.mask)[0]
  def random_choice(self):
    if self.mask is None:
      raise ValueError('Selection object is not initialized!')
    choice = np.random.choice(self.indices)
    # OutData['val'] = self.select_from[choice]
    # OutData['index'] = choice
    return choice



