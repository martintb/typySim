import numpy as np
from typySim.molecule import DummyMolecule
from typySim.core import BondList
from  typySim.core import PairTable

class System(object):
  ''' 
  The top-level container/controller for all simulation data. Nearly all
  worker objects contain pointers back to this object in order to access
  the beads coordinates, types, etc. The user is responsible for "building"
  the system, adding a simulation :class:`Box`, and loading an appropriate 
  :class:`Engine`.


  Attributes
  ----------
  nbeads : int
      Current number of beads in the simulation

  x,y,z : float ndarrary, size (nbeads)
      Arrays of bead cartesian coordinate positions

  types : int ndarray, size (nbeads)
      Array of bead types. These types are used to identifying interaction
      types, and, most importantly, how the Monte Carlo moves are applied to
      a molecule. They should be integer numbers starting from 0. 

  molecules : list, variable size
      A list containing all molecules that have been added to the system. 

  molecule_map : object ndarray, size (nbeads)
      A dictionary mapping between bead index and the molecule that is 
      associated with that bead index. A consequence of this data_structure 
      is that each bead index can only belong to one molecule at a time.

  bonds :  int ndarray, size (nbeads,max_nbonds)
      A list of bonds  for each beads. Each bond list should be sorted with
      -1 values signifying that no bond is present. 

  box : object
      A :class:`Box` object which handles all periodic wrapping and 
      distance calculation.

  engine : object
      A :class:`Engine` object which handles updating and modifying the system a
      according to a predefined scheme. 

  DummyMolecule : object
      An instance of the :class:`DummyMolecule`. The intention is that this
      object be used as a singleton, sentinel for testing purposed. This means
      that it be directly tested against using "is" constructs. 

  '''
  def __init__(self):
    #Basic containers and counters
    self.nbeads         = 0
    self.x              = np.array([],dtype=np.float)
    self.y              = np.array([],dtype=np.float)
    self.z              = np.array([],dtype=np.float)
    self.types          = np.array([],dtype=np.int)
    self.molecule_map   = np.array([],dtype=object)
    self.molecules      = list()
    self.computes       = list()
    self.bonds = BondList()

    self.trial_x      = None
    self.trial_y      = None
    self.trial_z      = None
    self.trial_types  = None
    self.trial_bond_pairlist  = None
    self.nbeads_trial = 0

    # Placeholders for worker objects to be set later. These variables have underscores
    # as their setting/getting behavior will be handled via the @property construct in 
    # python. This will allow more control over how the worker objects are initialized,
    # without having load_box and load_engine methods. 
    self._box = None 
    self._engine = None
    self.neighbor_list = None


    #Placeholders for Table objects with contain the parameters associated with bonded and
    # nonbonded interactions
    self.NonBondedTable = None
    self.BondedTable = None
    
    # sentinel dummy molecule object. A reference to this molecule will be placed 
    # in the molecules list when that bead is not "owned" by a molecule
    self.DummyMolecule = DummyMolecule()

    self.max_nbonds = 5 #arbitrary maximum on the number of bonds a single atom can have
  def set_trial_move(self,x,y,z,types,bonds):
    try:
      nbeads = len(x)
    except TypeError:
      nbeads = 1
      x = [x]
      y = [y]
      z = [z]
      types = [types]

    if self.nbeads_trial == nbeads:
      np.copyto(self.trial_x,x)
      np.copyto(self.trial_y,y)
      np.copyto(self.trial_z,z)
      np.copyto(self.trial_types,types)
      self.trial_bond_pairlist = np.array(bonds)
    else:
      self.nbeads_trial = nbeads
      self.trial_x = np.array(x)
      self.trial_y = np.array(y)
      self.trial_z = np.array(z)
      self.trial_types = np.array(types)
      self.trial_bond_pairlist = np.array(bonds)
  def append_trial_move(self):
    return self.add_beads(x=self.trial_x,
                          y=self.trial_y,
                          z=self.trial_z,
                          types=self.trial_types,
                          bonds=self.trial_bond_pairlist)
  @property
  def positions(self):
    return np.array([self.x,self.y,self.z],dtype=np.float).T
  @positions.setter
  def positions(self,pos):
    self.x = pos[:,0]
    self.y = pos[:,1]
    self.z = pos[:,2]
  @property
  def engine(self):
    return self._engine
  @engine.setter
  def engine(self,engine):
    if self._engine is not None:
      warngings.warn('--> Replacing currently loaded engine: {}'.format(self._engine.name))
    engine.system = self
    self._engine = engine
  @property
  def box(self):
    return self._box
  @box.setter
  def box(self,box):
    box.system = self
    self._box = box
    if self._box.neighbor_list is not None:
      self.neighbor_list = self._box.neighbor_list
  def reset_all_molecules(self,index_mapping=None):
    ''' Calls the :func:`reset` function on all molecules in 
        :attr:`self.molecules` so that their position/type masked arrays can be updated. 


        Parameters
        ----------
        index_mapping : dict, *Optional*
            If provided, the mapping dictionary will be passed to each molecules :func:`reset`
            function. The molecules can use this dictionary to update their indices.
        
    '''
    remove_list = []
    for i,mol in enumerate(self.molecules):
      keep_molecule = mol.reset(index_mapping)
      if not keep_molecule:
        remove_list.append(mol)

    for mol in remove_list:
      self.remove_molecule(mol,keep_beads=False)
  def add_beads(self, x,y,z,types, bonds=None,bond_shift=False):
    '''
    Primary method for adding new beads to the :class:`System`. Note that this method does
    not modify the molecule list or the molecules contained therein!

    Parameters
    ----------
    x,y,z : float ndarray, size (N), *Required*
        Array of positions to be added to system. No overlap or sanity check
        is carried out of these positions. The positions are also not 
        immediately wrapped, but this is likely to happen during :class:`Engine` 
        execution.

    types: float ndarray, size (N), *Required*
        Array of types to be added to the system.  

    bonds: list of pairs, *Optional*
        A list of bead pairs which are connected by a bond. Each bond will be added in
        forward and reverse direction (i-bonded-j and j-bonded-i). The :class:`System` class 
        uses sets to store bond pairs so duplicate bonds are not a worry.

    '''
    nx = len(x)
    ny = len(y)
    nz = len(z)
    nt = len(types)
    if not (nx==ny==nz==nt):
      raise ValueError(
                    '''Coordinate and types arrays do not have matching sizes!
                       num_x:{}
                       num_y:{}
                       num_z:{}
                       num_types:{}'''.format(nx,ny,nz,nt)
                      )
    new_nbeads = nx
    self.x     = np.append(self.x,x)
    self.y     = np.append(self.y,y)
    self.z     = np.append(self.z,z)
    self.types = np.append(self.types,types)
    
    #the molecule_map and bonds data structures must be initialized regardless of whether the beads
    #have bonds or are associated with molecules. 
    self.molecule_map = np.append(self.molecule_map,[self.DummyMolecule for i in range(new_nbeads)])

    self.bonds.expand(new_nbeads)
    if bonds is not None:
      if bond_shift == True:
        bond_shift_value = self.nbeads
      else:
        bond_shift_value = 0
      for i,j in bonds:
        self.bonds.add(i,j,bond_shift_value)

    #update the neighbor_list
    if self.neighbor_list is not None:
      for ix,iy,iz in zip(x,y,z):
        self.neighbor_list.insert_bead(-1,ix,iy,iz)

    #We'll return a list of the new indices that we've just added to the system
    #This way the user can immediately act on these added beads e.g. set their molecule
    new_indices = range(self.nbeads,self.nbeads+new_nbeads)
    self.nbeads += new_nbeads
    return new_indices
  # @profile
  def remove_beads(self, indices):
    '''
    Primary method for removing beads to the :class:`System`. Note that this method
    does not alter the molecule list or molecules contained therin!

    Parameters
    ----------
    indices : int list, size (N), *Required*
        List of indices to be removed from the System.
    '''

    #create index map. These will be needed for
    #converting molecules indices

    old2new = np.full((self.nbeads),-1,dtype=np.int)
    new2old = np.full((self.nbeads-len(indices)),-1,dtype=np.int)
    newi = 0
    for oldi in range(self.nbeads):
      if oldi in indices:
        old2new[oldi] = -1
      else:
        old2new[oldi] = newi
        new2old[newi] = oldi
        newi+=1

    removed_nbeads = len(indices)
    self.nbeads -= removed_nbeads

    self.x            = np.delete(self.x,indices)
    self.y            = np.delete(self.y,indices)
    self.z            = np.delete(self.z,indices)
    self.types        = np.delete(self.types,indices)
    self.molecule_map = np.delete(self.molecule_map,indices)

    self.bonds.shrink(indices)
    self.bonds.remap(old2new)

    #update the neighbor_list
    if self.neighbor_list is not None:
      # If the neighborlist was resized as we were removing the beads from the cells
      # we would need to conditionally map the indices as we were removing beads. By
      # making it a two step process we remove the beads from the linked-list chains,
      # and then go about resizing the list so that it matches the system. 
      for ix in indices:
        self.neighbor_list.remove_bead(ix)
      self.neighbor_list.shrink(indices)

    return {'old2new':old2new,'new2old':new2old,'n2o':new2old,'o2n':old2new}
  def add_molecule(self,molecule,**kwargs):
    ''' 
    Adds a molecule reference object to the system. If kwargs are supplied, new beads are
    created and the system containers are resized. This method handles the two possible 
    cases with adding a molecule to the system:

    Case 1: The molecule is made of up new beads and the system containers (x,y,z,types,etc)
    need to be appropriately expanded and updated. In this case, a called to 
    :func:`self.add_beads` needs to be made. 
    This situation is activiated by passing the appropriate keyword arguments 
    to this function which are relayed to :func:`self.add_beads`

    Case 2: The molecule is made up of beads already added to the system. In this case, the 
    molecule simply needs to be updated with a reference to this system and added to
    this systems molecule list. 
    '''
    if kwargs:
      molecule.indices = self.add_beads(**kwargs)
    molecule.system = self
    self.molecule_map[molecule.indices] = molecule
    self.molecules.append(molecule)

    # # Since new beads have been (possibly) added, it's necessary 
    # # to recreate all of the masked arrays of all of the molecules
    # if kwargs:
    #   self.reset_all_molecules()
    # else:
    self.molecules[-1].reset()
  def remove_molecule(self,molecule=None,index=None,remove_beads=False):
    ''' 
    Removes a reference molecule, and possible the beads that it references,
    from the system.

    Parameters
    ----------
    molecule : object, *optional*
        Molecule object to be removed. If supplied, the molecule list is searched until a
        molecule object with a matching reference (i.e. a matching :func:`id`) is found. This 
        found object is then removed.

    index : int, *optional*
        Index of molecule in :attr:`self.molecules` to be removed.
      
    remove_beads : bool, *optional*
        If True, :func:`self.remove_beads` is called on all indices inside the molecule being
        removed. If false, only the molecule reference object is removed. 

    '''
    if (molecule is None) and (index is None):
      raise ValueError('Must supply either index or molecule object to be removed!')
    if (molecule is not None) and (index is not None):
      raise ValueError('Both an index and a molecule object is supplied. Unsure of which to remove.')
    
    if molecule is not None:
      for i,mol in enumerate(self.molecules):
        if mol is molecule:
          index = i
          break

    if index is None:
      raise ValueError('Supplied molecule was not found in the molecule list.')

    removed_mol = self.molecules.pop(index)


    if remove_beads is True:
      index_mapping = self.remove_beads(removed_mol.indices)
      self.reset_all_molecules(index_mapping['old2new'])
    else:
      for i in range(self.nbeads):
        if self.molecule_map[i] is removed_mol:
          self.molecule_map[i] = self.DummyMolecule
    return removed_mol
  def get_compute(self,name):
    for compute in self.computes:
      if compute.name==name:
        return compute
    raise ValueError('Could not find a compute with name {}'.format(name))
  def add_compute(self,compute):
    self.computes.append(compute)


    

