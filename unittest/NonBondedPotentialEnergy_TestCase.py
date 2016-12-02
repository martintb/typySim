import unittest
import numpy as np
from math import sqrt
from typySim.compute import NonBondedPotentialEnergy
from typySim.core import System,Box,PairTable
from itertools import cycle
import pdb; ist = pdb.set_trace


class NonBondedPotentialEnergy_TestCase(unittest.TestCase):
  def hard_sphere(self,r,epsilon,sigma,rcut):
    if r<sigma:
      return 1e9
    else:
      return 0
  def lennard_jones(self,r,epsilon,sigma,rcut):
    if r<rcut:
      return 4*epsilon*((sigma/r)**(12.0) - (sigma/r)**(6.0))
    else:
      return 0
  def calc_all_potential(self,x,y,z,t,lx,ly,lz,PT):
    N = len(x)
    U = 0
    for i in range(N):
      for j in range(N):
        if (j<=i):
          continue
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
        if dist==0:
          print 'i',i,x[i],y[i],z[i],t[i],'j',j,x[j],y[j],z[j],t[j]
        eps = PT['epsilon',t[i],t[j]]
        sig = PT['sigma',t[i],t[j]]
        rcut = PT['rcut',t[i],t[j]]
        if PT['potential',t[i],t[j]]=='HardSphere':
          U+=self.hard_sphere(dist,eps,sig,rcut)
        elif PT['potential',t[i],t[j]]=='LennardJones':
          U+=self.lennard_jones(dist,eps,sig,rcut)
        else:
          raise ValueError('Not set up for this potential!')
    return U
  def make_rubix(self,N=1,dx=1,types=None):
    if types is None:
      types = cycle([0])
    else:
      types = cycle(types)
    #Basic sanity check that beads are being added
    x = []; y = []; z = []; t = []
    for i in np.arange(-N,N,dx):
      for j in np.arange(-N,N,dx):
        for k in np.arange(-N,N,dx):
          x.append(i)
          y.append(j)
          z.append(k)
          t.append(types.next())
    return x,y,z,t
  def base_test(self,PT,lx,ly,lz,types=None,N=3,box_resize=None):
    x,y,z,t = self.make_rubix(N=N,types=types)
    U0 = self.calc_all_potential(x,y,z,t,lx,ly,lz,PT)
    

    system = System()
    system.add_beads(x=x,y=y,z=z,types=t)
    system.box = Box(cell_grid=(3,3,3))
    system.box.lx = lx
    system.box.ly = ly
    system.box.lz = lz

    xarray = np.array(x,dtype=np.float)
    yarray = np.array(y,dtype=np.float)
    zarray = np.array(z,dtype=np.float)
    system.box.neighbor_list.build_nlist(xarray,yarray,zarray,True)

    system.NonBondedTable = PT
    PE = NonBondedPotentialEnergy(system)
    U1 = PE.compute()
    U2 = PE.compute(ignore_neighbor_list=True)
    if box_resize is None:
      return U0,U1,U2
    system.box.lx = box_resize[0]
    system.box.ly = box_resize[1]
    system.box.lz = box_resize[2]
    U3 = self.calc_all_potential(x,y,z,t,system.box.lx,system.box.ly,system.box.lz,PT)
    U4 = PE.compute()
    U5 = PE.compute(ignore_neighbor_list=True)
    return U0,U1,U2,U3,U4,U5
  def test_nopbc_homog_HS_sigma10(self):
    lx = 100
    ly = 100
    lz = 100

    PT = PairTable(types=['P1'],parms=['epsilon','rcut','sigma','potential'])
    PT.setUnsetValues('epsilon',1.0)
    PT.setUnsetValues('sigma',1.0)
    PT.setUnsetValues('rcut',2.5)
    PT.setUnsetValues('potential','HardSphere')

    U0,U1,U2 = self.base_test(PT,lx,ly,lz)
    self.assertAlmostEqual(U0,U1)
    self.assertAlmostEqual(U0,U2)
  def test_nopbc_homog_LJ_sigma10(self):
    lx = 100
    ly = 100
    lz = 100

    PT = PairTable(types=['P1'],parms=['epsilon','rcut','sigma','potential'])
    PT.setUnsetValues('epsilon',1.0)
    PT.setUnsetValues('sigma',1.0)
    PT.setUnsetValues('rcut',2.5)
    PT.setUnsetValues('potential','LennardJones')

    U0,U1,U2 = self.base_test(PT,lx,ly,lz)
    self.assertAlmostEqual(U0,U1)
    self.assertAlmostEqual(U0,U2)
  def test_nopbc_homog_HS_sigma155(self):
    lx = 100
    ly = 100
    lz = 100

    PT = PairTable(types=['P1'],parms=['epsilon','rcut','sigma','potential'])
    PT.setUnsetValues('epsilon',1.0)
    PT.setUnsetValues('sigma',1.55)
    PT.setUnsetValues('rcut',2.5)
    PT.setUnsetValues('potential','HardSphere')

    U0,U1,U2 = self.base_test(PT,lx,ly,lz)
    self.assertAlmostEqual(U0,U1)
    self.assertAlmostEqual(U0,U2)
  def test_pbc_homog_HS_sigma155(self):
    lx = 10
    ly = 11
    lz = 10

    PT = PairTable(types=['P1'],parms=['epsilon','rcut','sigma','potential'])
    PT.setUnsetValues('epsilon',1.0)
    PT.setUnsetValues('sigma',1.55)
    PT.setUnsetValues('rcut',2.5)
    PT.setUnsetValues('potential','HardSphere')

    U0,U1,U2 = self.base_test(PT,lx,ly,lz)
    self.assertAlmostEqual(U0,U1)
    self.assertAlmostEqual(U0,U2)
  def test_pbc_homog_LJ_sigma10(self):
    lx = 10
    ly = 11
    lz = 10

    PT = PairTable(types=['P1'],parms=['epsilon','rcut','sigma','potential'])
    PT.setUnsetValues('epsilon',1.0)
    PT.setUnsetValues('sigma',1.0)
    PT.setUnsetValues('rcut',2.5)
    PT.setUnsetValues('potential','LennardJones')

    U0,U1,U2 = self.base_test(PT,lx,ly,lz)
    self.assertAlmostEqual(U0,U1)
    self.assertAlmostEqual(U0,U2)
  def test_pbc_mixed_LJ(self):
    lx = 8
    ly = 9
    lz = 8

    PT = PairTable(types=['P1','P2'],parms=['epsilon','rcut','sigma','potential'])
    PT['epsilon','P1','P1'] = 2.0
    PT['epsilon','P1','P2'] = 1.0
    PT['epsilon','P2','P2'] = 1.25
    PT['sigma','P1','P1'] = 1.0
    PT['sigma','P1','P2'] = 1.1
    PT['sigma','P2','P2'] = 1.25
    PT['potential','P1','P2'] = 'HardSphere'
    PT.setUnsetValues('potential','LennardJones')
    PT.setUnsetValues('rcut',2.5)

    U0,U1,U2,U3,U4,U5 = self.base_test(PT,lx,ly,lz,N=4,box_resize=(100,100,100))
    self.assertAlmostEqual(U0,U1)
    self.assertAlmostEqual(U0,U2)
    self.assertAlmostEqual(U3,U4)
    self.assertAlmostEqual(U3,U5)
    self.assertNotAlmostEqual(U0,U3)




if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(NonBondedPotentialEnergy_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


