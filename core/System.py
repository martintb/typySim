import numpy as np
from typySim.molecule import DummyMolecule
from  PairTable import PairTable

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

  bonds :  list, size (nbeads)
      A list of bonds  for each beads. The list datatype is used
      because the number of bonds for each beads can vary.

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
    self.bonds          = np.array([],dtype=object)
    self.molecules      = list()
    #self. compute

    # Placeholders for worker objects to be set later. These variables have underscores
    # as their setting/getting behavior will be handled via the @property construct in 
    # python. This will allow more control over how the worker objects are initialized,
    # without having load_box and load_engine methods. 
    self._box = None 
    self._engine = None
    self.neighbor_list = None
    
    # sentinel dummy molecule object. A reference to this molecule will be placed 
    # in the molecules list when that bead is not "owned" by a molecule
    self.DummyMolecule = DummyMolecule()
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
    for mol in self.molecules:
      mol.reset(index_mapping)
  def add_beads(self, x,y,z,types, bonds=None):
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
    self.bonds        = np.append(self.bonds,       [set()              for i in range(new_nbeads)])
    if bonds is not None:
      for i,j in bonds:
        self.add_bond(i,j)

    #We'll return a list of the new indices that we've just added to the system
    #This way the user can immediately act on these added beads e.g. set their molecule
    new_indices = range(self.nbeads,self.nbeads+new_nbeads)
    self.nbeads += new_nbeads
    return new_indices
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
    old2new = {}
    new2old = {}
    newi = 0
    for oldi in range(self.nbeads):
      if oldi in indices:
        old2new[oldi] = None
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
    self.bonds        = np.delete(self.bonds,indices)

    # Processing the bonds is a bit of a pain. We need to remove all
    # references to removed indices while also mapping from old to 
    # new index numbers.
    new_bonds = []
    for old_bonded_list in self.bonds:
      new_bonded_list = []
      for i in old_bonded_list:
        if old2new[i] is not None:
          new_bonded_list.append(old2new[i])
      new_bonds.append(set(new_bonded_list))
    self.bonds = np.array(new_bonds,dtype=object)

    return {'old2new':old2new,'new2old':new2old,'n2o':new2old,'o2n':old2new}
  def add_bond(self,i,j,shift=None):
    '''
    Parameters
    ----------
    i,j : int, *required*
        Bond indices to be added to the :class:`System`.

        Bond indices passed as :class:`int` are assumed to be relative to the group of beads 
        being added.This means the passing [[0,1]] will add a bond between the first and 
        second beads **of the beads currently being added** and not the first and second 
        beads of the system. This effect is achieved by shifting the passed bond indices 
        by self.nbeads.

        Bond indices passed as :class:`str` indicate that you wish to bond against an index
        that has already been added o the system. These bond indices are made positive, but 
        left unshifted. 

        .. Warning::
          No sanity checks are carried out on the bond indices. Make sure you don't attempt
          to bond an bead index that doesn't exist in the system!

        .. Warning::
          This approach is hacky and bad and alternative ideas would be accepted. The basic problem
          is that a passed molecule needs to be able to specify all new internal bonds and any bonds
          to existing beads simultaneously. Bonds between new and old beads need to be possible. 
    shift : bool, *optional*
        Value to use as the shift-value if the beads are passed as :class:`int`. The default value
        is :func:`self.nbeads`.
          
          
     '''
    if shift is None:
      shift = self.nbeads

    if isinstance(i,basestring):
      bondi = int(i)
    else:
      bondi = i+shift

    if isinstance(j,basestring):
      bondj = int(j)
    else:
      bondj = j+shift

    self.bonds[bondi].add(bondj)
    self.bonds[bondj].add(bondi)
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

    # Since new beads have been (possibly) added, it's necessary 
    # to recreate all of the masked arrays of all of the molecules
    self.reset_all_molecules()
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
        if id(mol)==id(molecule):
          index = i
          break

    if index is None:
      raise ValueError('Supplied molecule was not found in the molecule list.')

    removed_mol = self.molecules.pop(index)

    if remove_beads is True:
      index_mapping = self.remove_beads(removed_mol.indices)
      self.reset_all_molecules(index_mapping['old2new'])
    return removed_mol
    

