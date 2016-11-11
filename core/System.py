import numpy as np

class System(object):
  ''' 
  safe2raw:
  ---------
  The indices array is essentially a map from "safe_indices" to "raw_indices". This array can never shrink in size and when an atom is removed from the system, the value at that atoms index becomes None (or zero...). In essence, each element in the indices array is persistent for the duration of the simulation. 

  safe_index: 
  ----------
  A safe_index refers to an element in the indices array. As these elements can never be deleted (only reassigned), these indices should refer to the same "atom" for the lifetime of the system object. 

  raw index:  
  ----------
  A raw_index refers to an element in almost all other arrays. After atom addition, deletion, or resorting, the value of these indices may no longer have meaning. 
  
  '''
  def __init__(self):
    '''
    The first bead is a dummy atom that is used as a placeholder for removed beads
    '''
    self.natoms         = 0
    self.safe2raw       = np.array([0])
    self.positions      = np.array([[0.123,0.456,0.789]],dtype=float)
    self.types          = np.array([-1],dtype=np.int)
    self.bonds          = list()
    self.molecules      = list()
    self.next_index     = 1
  def reset_all(self):
    for mol in self.molecules:
      mol.reset()
  def add_atoms(self,**kwargs):
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

    if 'types' in kwargs:
      self.types = np.append(self.types,kwargs['types'])

    if 'bonds' in kwargs:
      self.bonds = np.append(self.bonds,kwargs['bonds'])

    return new_indices
  def add_molecule(self,molecule,**mol_kwargs):
    molData = molecule.build(**mol_kwargs)
    molecule.system = self
    molecule.indices = self.add_atoms(**molData)
    self.molecules.append(molecule)
    

