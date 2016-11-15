import warnings
import numpy as np
from Box import Box
from ..molecule import DummyMolecule

class System(object):
  ''' 
  The top-level container/controller for all simulation data

  Attributes
  ----------
  :ivar nbeads: current  number of beads in the simulation

  :ivar positions: array of bead x,y,z positions

  :ivar types: array of integer bead types

  :ivar molecules: list of current molecule objects

  :ivar molecule_map: mapping between bead index and primary molecule
                      a consequence of this data_structure is that each 
                      bead can only belong to one molecule at a time.

  :ivar bonds: list of bonds  for each atom. The list datatype is used
               because the number of bonds for each atom can vary.

  :ivar box: box object which handles all periodic wrapping and distance
             calculation.

  :ivar DummyMolecule: used as a placeholder "sentinel" that can
                       be directly tested against using "is" constructs as shown
                       below. No other instances of thi
                       >>> system.molecules[0] is system.DummyMolecule # True/False
  '''
  def __init__(self,cell_grid=None):
    self.nbeads         = 0
    self.positions      = None
    self.types          = np.array([],dtype=np.int)
    self.molecules      = list()
    self.molecule_map   = np.array([])
    self.bonds          = list()

    self.box = Box(system=self,cell_grid=cell_grid)
    self.computes = list()
    self.engine = None

    self.DummyMolecule = DummyMolecule()
  def reset_all_molecules(self):
    '''
    Attempts to reset all molecules so position/type masked arrays can be updated
    '''
    for mol in self.molecules:
      mol.reset()
  def add_beads(self,**kwargs):
    '''
    Primary method for adding new beads to the System.

    Parameters
    ----------
    :param positions: `Required`
    :param types: `Required` 
    '''
    try:
      new_positions = kwargs['positions']
    except KeyError:
      raise ValueError('System.add_beads called without positions kwarg!')
    if self.positions is None:
      self.positions = np.array(new_positions,dtype=float)
    else:
      try:
        self.positions = np.append(self.positions,new_positions,axis=0)
      except ValueError:
        import pdb; pdb.set_trace()


    try:
      new_types = kwargs['types']
    except KeyError:
      raise ValueError('System.add_beads called without types kwarg!')
    self.types = np.append(self.types,new_types)

    nbeads_new = len(kwargs['positions'])
    new_indices = range(self.nbeads,self.nbeads+nbeads_new)
    self.bonds.extend([set() for i in range(nbeads_new)])
    self.molecule_map = np.append(self.molecule_map,[self.DummyMolecule for i in range(nbeads_new)])

    if 'bonds' in kwargs:
      for i,j in kwargs['bonds']:
        self.add_bond(i,j)

    self.nbeads += nbeads_new
    return set(new_indices)
  def add_bond(self,i,j):
    self.bonds[i].add(j)
    self.bonds[j].add(i)
  def add_molecule(self,molecule,molData=None):
    if molData is not None:
      molecule.indices = self.add_beads(**molData)
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
      warngings.warn('--> Replacing currently loaded engine: {}'.format(self.engine.name))
    engine.system = self
    self.engine = engine
    

