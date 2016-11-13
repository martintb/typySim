from Compute import Compute


class PotentialEnergy(Compute):
  def __init__(self):
    super(PotentialEnergy,self).__init__() #call parent classes constructor
    self.name='PotentialEnergy'
  def compute(self):
    if self.system is None:
      raise ValueError('{} is not associated with a system! Cannot compute!'.format(self.name))

