import numpy as np
from Box import Box

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
    self.bonds          = np.array([0,0,0,0],dtype=np.int) #each bead can have up to four bonds, zero means no bond
    self.molecules      = list()
    self.next_index     = 1
    self.box = Box()
    self.computes = list()
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

    if 'bonds' in kwargs:
      self.bonds = np.append(self.bonds,kwargs['bonds'])

    return new_indices
  def add_molecule(self,molecule,build=True,**kwargs):
    if build==True:
      molData = molecule.build(**kwargs)
    elif build==False:
      molData = kwargs
    molecule.system = self
    molecule.indices = self.add_atoms(**molData)
    self.molecules.append(molecule)
  def get_compute(self,name):
    for k,v in self.computes.items():
      if k==name:
        return v
    

