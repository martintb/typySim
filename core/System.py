import numpy as np
from Box import Box
from ..molecule import DummyMolecule

class System(object):
  def __init__(self):
    '''
    The first bead is a dummy atom that is used as a placeholder for removed beads
    '''
    self.DummyMolecule = DummyMolecule() #sentinel DummyMolecule
    self.natoms         = 0
    self.positions      = None
    self.types          = np.array([],dtype=np.int)
    self.molecules      = list()
    self.molecule_map   = np.array([])
    self.bonds          = [] #dummy atom has no bonds

    self.box = Box(system=self)
    self.computes = list()
    self.engine = None
  def reset_all_molecules(self):
    for mol in self.molecules:
      mol.reset()
  def add_atoms(self,**kwargs):
    try:
      if self.positions is None:
        self.positions = np.array(kwargs['positions'],dtype=float)
      else:
        self.positions = np.array(np.append(self.positions,kwargs['positions'],axis=0))
    except KeyError:
      raise ValueError('System.add_atoms called without positions kwarg!')


    try:
      self.types = np.array(np.append(self.types,kwargs['types']))
    except KeyError:
      raise ValueError('System.add_atoms called without types kwarg!')

    natoms_new = len(kwargs['positions'])
    new_indices = range(self.natoms,self.natoms+natoms_new)
    self.bonds.extend([set() for i in range(natoms_new)])
    self.molecule_map = np.append(self.molecule_map,[self.DummyMolecule for i in range(natoms_new)])

    if 'bonds' in kwargs:
      for bondi,bondj in kwargs['bonds']:
        # if bondi>=0:
        #   bondi+=self.natoms
        # if bondj>=0:
        #   bondj+=self.natoms
        # bondi = abs(bondi) #handle the negative bonds
        # bondj = abs(bondj) #handle the negative bonds
        self.bonds[bondi].add(bondj)
        self.bonds[bondj].add(bondi)
        
    self.natoms += natoms_new
    return set(new_indices)
  def add_molecule(self,molecule,molData=None):
    if molData is not None:
      molecule.indices = self.add_atoms(**molData)
    molecule.system = self
    self.molecule_map[list(molecule.indices)] = molecule
    self.molecules.append(molecule)
    self.reset_all_molecules()
  def get_compute(self,name):
    for k,v in self.computes.items():
      if k==name:
        return v
  def add_engine(self,engine):
    if self.engine is not None:
      raise Warning('Replacing currently loaded engine: {}'.format(self.engine.name))
    engine.system = self
    self.engine = engine
    

