import numpy as np
import logging

class MoleculeCalculator(object):
  '''Compositional class that holds all "calculator" methods

  This class is intended to be instantiated only by the molecule base class
  and then used as an interface by all sub-molecules for all common calculations.

  '''
  def __init__(self,molecule):
    self.molecule = molecule
  def center_of_mass(self):
    '''
    .. Warning::
      This class returns the unwrapped center of mass. For free floating molecules,
      this position could be many box lengths away from the central image. The current
      :class:`typySim.core.Box` implementation does not handle long range unwrapping!

    '''

    ux = self.molecule.x.compressed()
    uy = self.molecule.y.compressed()
    uz = self.molecule.z.compressed()

    box = self.molecule.system.box
    if box is not None:
      #Unwrap using image flags
      imx = self.molecule.imx.compressed()
      imy = self.molecule.imy.compressed()
      imz = self.molecule.imz.compressed()

      factor_x = imx*box.lx
      factor_y = imy*box.ly
      factor_z = imz*box.lz
    else:
      factor_x = 0
      factor_y = 0
      factor_z = 0

    com_x = np.mean(ux - factor_x)
    com_y = np.mean(uy - factor_y)
    com_z = np.mean(uz - factor_z)
    com = np.array([com_x,com_y,com_z])

    self.molecule.properties['center_of_mass'] = com
    return com




