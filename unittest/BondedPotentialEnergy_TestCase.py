import unittest
import numpy as np
from math import sqrt
from math import log
from typySim.compute import BondedPotentialEnergy
from typySim.core import System,Box,PairTable
from itertools import cycle
import ipdb as pdb; ist = pdb.set_trace

import struct
def test(num,n=8):
  return ''.join(bin(ord(c)).replace('0b', '').rjust(n, '0') for c in struct.pack('!f', num))


class BondedPotentialEnergy_TestCase(unittest.TestCase):
  def harmonic(self,r,k,r0):
    return 0.5 * k * (r-r0)**(2.0)
  def FENE(self,r,k,r0):
    return 0.5*k*r0*r0*log(1-(r/r0)*(r/r0))
  def calc_all_potential(self,bonds,x,y,z,t,lx,ly,lz,PT):
    U = 0
    
    for i,j in bonds:
      dx = abs(x[i] - x[j])
      dy = abs(y[i] - y[j])
      dz = abs(z[i] - z[j])
      if dx>lx/2.0:
        dx-=lx
      if dy>ly/2.0:
        dy-=ly
      if dz>lz/2.0:
        dz-=lz
      dist = sqrt(dx*dx+dy*dy+dz*dz)
      k = PT['k',t[i],t[j]]
      r0 = PT['r0',t[i],t[j]]
      if PT['potential',t[i],t[j]]=='Harmonic':
        U+=self.harmonic(dist,k,r0)
      elif PT['potential',t[i],t[j]]=='FENE':
        U+=self.FENE(dist,k,r0)
      else:
        raise ValueError('Not set up for this potential!')
    return U
  def base_test(self,lx,ly,lz,N,dz,potential,r0):
    x = [0.0]*N
    y = [0.0]*N
    z = [dz*i-lz/2.0 for i in range(N)]
    t = [0]*N
    bonds = [[i,j] for i,j in zip(range(N-1),range(1,N))]


    NBPT = PairTable(types=['P1'],parms=['epsilon','rcut','sigma','potential'])
    NBPT.setUnsetValues('epsilon',1.0)
    NBPT.setUnsetValues('sigma',1.0)
    NBPT.setUnsetValues('rcut',2.5)
    NBPT.setUnsetValues('potential','HardSphere')

    BPT = PairTable(types=['P1'],parms=['k','r0','potential'])
    BPT.setUnsetValues('k',1.0)
    BPT.setUnsetValues('r0',r0)
    BPT.setUnsetValues('potential',potential)

    system = System()
    system.add_beads(x=x,y=y,z=z,types=t,bonds=bonds)
    system.box = Box(cell_grid=(3,3,3))
    system.box.lx = lx
    system.box.ly = ly
    system.box.lz = lz
    system.BondedTable = BPT
    system.NonBondedTable = NBPT


    x,y,z = system.box.numpy_wrap_position(x=x,y=y,z=z).T
    U0 = self.calc_all_potential(bonds,x,y,z,t,lx,ly,lz,BPT)

    xarray = np.array(x,dtype=np.float)
    yarray = np.array(y,dtype=np.float)
    zarray = np.array(z,dtype=np.float)
    system.box.neighbor_list.build_nlist(xarray,yarray,zarray,True)

    PE = BondedPotentialEnergy(system)
    U1 = PE.compute()

    return U0,U1
  def test_harmonic(self):
    lx = 8
    ly = 11
    lz = 9
    N = 10
    dz = 1.1
    r0 = 1.0
    potential = 'Harmonic'
    U0,U1 = self.base_test(lx,ly,lz,N,dz,potential,r0)
    self.assertAlmostEqual(U0,U1,delta=0.001)
  def test_FENE(self):
    lx = 8
    ly = 11
    lz = 9
    N = 10
    dz = 1.1
    r0 = 2.0
    potential = 'FENE'
    U0,U1 = self.base_test(lx,ly,lz,N,dz,potential,r0)
    self.assertAlmostEqual(U0,U1,delta=0.001)



if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(BondedPotentialEnergy_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


