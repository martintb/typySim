import numpy as np
from Box import Box
from ..molecule import DummyMolecule

class System(object):
  ''' 
  safe2raw:
  ---------
  The indices array is essentially a map from "safe_indices" to "raw_indices". This array can never shrink in size and when an atom is removed from the system, the value at that atoms index becomes None (or zero...). In essence, each element in the indices array is persistent for the duration of the simulation. 
  
  '''
  def __init__(self):
    '''
    The first bead is a dummy atom that is used as a placeholder for removed beads
    '''
    self.natoms         = 0
    self.safe2raw       = np.array([0])
    self.positions      = np.array([[0.123,0.456,0.789]],dtype=np.float)
    self.types          = np.array([-1],dtype=np.int)
    self.molecules      = list()
    self.molecule_map   = np.array([DummyMolecule()])
    self.bonds          = [set()] #dummy atom has no bonds

    self.next_index     = 1
    self.box = Box()
    self.computes = list()
  def reset_all_molecules(self):
    for mol in self.molecules:
      mol.reset()
  def _add_atoms(self,**kwargs):
    try:
      self.positions = np.array(np.append(self.positions,kwargs['positions'],axis=0))
    except KeyError:
      raise ValueError('System.add_atoms called without positions kwarg!')

    try:
      self.types = np.array(np.append(self.types,kwargs['types']))
    except KeyError:
      raise ValueError('System.add_atoms called without types kwarg!')

    natoms_new = len(kwargs['positions'])
    new_indices = range(self.next_index,self.next_index+natoms_new)
    self.safe2raw  = np.append(self.safe2raw,new_indices)
    self.next_index = self.safe2raw[-1]+1
    self.natoms += natoms_new
    self.bonds.extend([set() for i in range(natoms_new)])
    self.molecule_map = np.append(self.molecule_map,[DummyMolecule() for i in range(natoms_new)])

    if 'bonds' in kwargs:
      for bondi,bondj in kwargs['bonds']:
        mapped_i = new_indices[bondi]
        mapped_j = new_indices[bondj]
        self.bonds[mapped_i].add(mapped_j)
        self.bonds[mapped_j].add(mapped_i)
        
    return new_indices
  def add_molecule(self,molecule,build=True,**kwargs):
    if build==True:
      molData = molecule.build(**kwargs)
    elif build==False:
      molData = kwargs
    molecule.system = self
    molecule.indices = self._add_atoms(**molData)
    self.molecule_map[molecule.indices] = molecule
    self.molecules.append(molecule)
  def get_compute(self,name):
    for k,v in self.computes.items():
      if k==name:
        return v
    

