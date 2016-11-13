from Molecule import Molecule

class Bead(Molecule):
  def __init__(self):
    super(Bead,self).__init__() #Need to call parent class' constructor
    self.name = 'Bead'
