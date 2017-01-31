import numpy as np
import logging
from typySim.core.cy.cyutil import n_closest
from typySim.core import Box

class MoleculeDistributor(object):
  '''
  '''
  def __init__(self,molecule):
    self.molecule = molecule
  def by_index(self,index_groups):
    '''
    Break this molecule up into smaller molecules of identical type.

    Parameters
    ----------
    index_groups : list, *Required*
        Index group is a two-dimensional list where len(index_groups) is 
        the number of sub-molecules to create and each list item is a list
        of indices to use for each molecule.

        .. Warning::
          All indices must be accounted for in the new molecules. An exception
          will be raised if this is not true.
    
    '''
    #Flatten list of lists into a 1-D set
    index_check = set([i for group in index_groups for i in group])
    if index_check!=set(self.molecule.indices):
      raise ValueError('Mismatch between passed indices and molecule indices.')

    new_molecules = []
    for group in index_groups:
      new_molecule = self.molecule.__class__()
      new_molecule.indices = group
      self.molecule.system.add_molecule(new_molecule)
      new_molecules.append(new_molecule)
    self.molecule.system.remove_molecule(molecule=self.molecule,remove_beads=False)
    return new_molecules
  def by_type_proximity(self,typeList):
    ''' Distribute a single molecule into multiple molecules based on proximity to a bead type
    '''
    if len(typeList)<2:
      raise ValueError('Need at least two types to distribute between!')

    box = self.molecule.system.box
    if box is None:
      box = Box(L=1000)

    x = self.molecule.x.compressed()
    y = self.molecule.y.compressed()
    z = self.molecule.z.compressed()
    types = self.molecule.types.compressed()
    indices = np.array(self.molecule.indices)


    # index_groups will be passed to self.distribute_by_index. It must be initialized with the
    # beads which are pre-grouped because they match a grouping type
    index_groups = []
    for i,typeVal in enumerate(typeList):
      index_groups.append(list(indices[types==typeVal]))

    # Need to combine all typeGroups. These beads are already grouped
    # as they define the groups themselves
    mask = (types==typeList[0])
    for typeVal in typeList[1:]:
      mask = np.logical_or(mask,types==typeVal)

    # Beads to be grouped
    x1Array  = x[~mask]
    y1Array  = y[~mask]
    z1Array  = z[~mask]
    indices1 = indices[~mask]

    # Beads which are pre-grouped because they match a proximity type
    x2Array = x[mask]
    y2Array = y[mask]
    z2Array = z[mask]
    types2  = types[mask]

    for i,(x1,y1,z1) in enumerate(zip(x1Array,y1Array,z1Array)):
      idex,dist = n_closest(1,x1,y1,z1,x2Array,y2Array,z2Array,box)
      for j,typeVal in enumerate(typeList):
        if types2[idex[0]] == typeVal:
          index_groups[j].append(indices1[i])
          break

    return self.by_index(index_groups)




