import unittest
import numpy as np
from typySim.core import System
from typySim import molecule


class Molecule_TestCase(unittest.TestCase):
  def test_add_single_molecule(self):
    ''' 
    Basic test to insure that a molecule is correctly added
    to a system
    '''

    x = range(8)
    y = range(8,16)
    z = range(16,24)
    types = [0,1,1,0]*2
    bonds = [[i,j] for i,j in zip(range(len(x)-1),range(1,len(x)))]
    chain = molecule.ChainSegment()

    system = System()
    system.add_molecule(chain,x=x,y=y,z=z,types=types,bonds=bonds)


    self.assertEqual(system.nbeads,len(x))
    self.assertIs(system.molecules[0],chain)
    self.assertIs(system.molecules[0].system,system)
    self.assertSetEqual(system.molecules[0].indices,set(range(8)))
    self.assertListEqual(list(system.x),x)
    self.assertListEqual(list(system.y),y)
    self.assertListEqual(list(system.z),z)
    self.assertListEqual(list(system.types),types)
    self.assertSetEqual(system.bonds[0],set([1]))
    self.assertSetEqual(system.bonds[1],set([0,2]))
    self.assertSetEqual(system.bonds[2],set([1,3]))
    self.assertSetEqual(system.bonds[3],set([2,4]))
    self.assertSetEqual(system.bonds[4],set([3,5]))
    self.assertSetEqual(system.bonds[5],set([4,6]))
    self.assertSetEqual(system.bonds[6],set([5,7]))
    self.assertSetEqual(system.bonds[7],set([6]))
  def test_remove_end_molecule(self):
    ''' 
    Basic test to insure that multiple molecules 
    can be added and removed to a system. Adds two ChainSegments
    to the system, checks to make sure that they were added correctly,
    and then removes the first chain
    '''

    x1 = range(4)
    y1 = range(4,8)
    z1 = range(8,12)
    types1 = [0,1,1,0]
    bonds1 = [[i,j] for i,j in zip(range(len(x1)-1),range(1,len(x1)))]
    chain1 = molecule.ChainSegment()

    x2 = range(10,14)
    y2 = range(14,18)
    z2 = range(18,22)
    types2= [2,3,3,2]
    bonds2 = [[i,j] for i,j in zip(range(len(x2)-1),range(1,len(x2)))]
    #add in bond between middle of first chain and beginning of second
    bonds2 += [['1',0]] 
    chain2 = molecule.ChainSegment()

    system = System()
    system.add_molecule(chain1,x=x1,y=y1,z=z1,types=types1,bonds=bonds1)
    system.add_molecule(chain2,x=x2,y=y2,z=z2,types=types2,bonds=bonds2)


    self.assertEqual(system.nbeads,len(x1+x2))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[1],chain2)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1].system,system)
    self.assertSetEqual(system.molecules[0].indices,set(range(4)))
    self.assertSetEqual(system.molecules[1].indices,set(range(4,8)))
    self.assertListEqual(list(system.x),x1+x2)
    self.assertListEqual(list(system.y),y1+y2)
    self.assertListEqual(list(system.z),z1+z2)
    self.assertListEqual(list(system.types),types1+types2)
    self.assertSetEqual(system.bonds[0],set([1]))
    self.assertSetEqual(system.bonds[1],set([0,2,4]))
    self.assertSetEqual(system.bonds[2],set([1,3]))
    self.assertSetEqual(system.bonds[3],set([2]))
    self.assertSetEqual(system.bonds[4],set([1,5]))
    self.assertSetEqual(system.bonds[5],set([4,6]))
    self.assertSetEqual(system.bonds[6],set([5,7]))
    self.assertSetEqual(system.bonds[7],set([6]))

    system.remove_molecule(molecule=chain1,remove_beads=True)
    self.assertEqual(system.nbeads,len(x2))
    self.assertIs(system.molecules[0],chain2)
    self.assertIs(system.molecules[0].system,system)
    self.assertSetEqual(system.molecules[0].indices,set(range(4)))
    self.assertListEqual(list(system.x),x2)
    self.assertListEqual(list(system.y),y2)
    self.assertListEqual(list(system.z),z2)
    self.assertListEqual(list(system.types),types2)
    self.assertSetEqual(system.bonds[0],set([1]))
    self.assertSetEqual(system.bonds[1],set([0,2]))
    self.assertSetEqual(system.bonds[2],set([1,3]))
    self.assertSetEqual(system.bonds[3],set([2]))
  def test_remove_middle_molecule(self):
    ''' 
    Basic test to insure that multiple molecules 
    can be added and removed to a system. Adds two ChainSegments
    to the system, checks to make sure that they were added correctly,
    and then removes the first chain
    '''

    x1 = range(4)
    y1 = range(4,8)
    z1 = range(8,12)
    types1 = [0,1,1,0]
    bonds1 = [[i,j] for i,j in zip(range(len(x1)-1),range(1,len(x1)))]
    chain1 = molecule.ChainSegment()

    x2 = range(10,14)
    y2 = range(14,18)
    z2 = range(18,22)
    types2= [2,3,3,2]
    bonds2 = [[i,j] for i,j in zip(range(len(x2)-1),range(1,len(x2)))]
    #add in bond between middle of first chain and beginning of second
    bonds2 += [['1',0]] 
    chain2 = molecule.ChainSegment()

    x3 = range(20,24)
    y3 = range(24,28)
    z3 = range(28,32)
    types3= [3,4,4,3]
    bonds3 = [[i,j] for i,j in zip(range(len(x3)-1),range(1,len(x3)))]
    #add in bond between middle of first chain and beginning of second
    chain3 = molecule.ChainSegment()

    system = System()
    system.add_molecule(chain1,x=x1,y=y1,z=z1,types=types1,bonds=bonds1)
    system.add_molecule(chain2,x=x2,y=y2,z=z2,types=types2,bonds=bonds2)
    system.add_molecule(chain3,x=x3,y=y3,z=z3,types=types3,bonds=bonds3)

    self.assertEqual(system.nbeads,len(x1+x2+x3))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[1],chain2)
    self.assertIs(system.molecules[2],chain3)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1].system,system)
    self.assertIs(system.molecules[2].system,system)
    self.assertSetEqual(system.molecules[0].indices,set(range(4)))
    self.assertSetEqual(system.molecules[1].indices,set(range(4,8)))
    self.assertSetEqual(system.molecules[2].indices,set(range(8,12)))
    self.assertListEqual(list(system.x),x1+x2+x3)
    self.assertListEqual(list(system.y),y1+y2+y3)
    self.assertListEqual(list(system.z),z1+z2+z3)
    self.assertListEqual(list(system.types),types1+types2+types3)
    self.assertSetEqual(system.bonds[0],set([1]))
    self.assertSetEqual(system.bonds[1],set([0,2,4]))
    self.assertSetEqual(system.bonds[2],set([1,3]))
    self.assertSetEqual(system.bonds[3],set([2]))
    self.assertSetEqual(system.bonds[4],set([1,5]))
    self.assertSetEqual(system.bonds[5],set([4,6]))
    self.assertSetEqual(system.bonds[6],set([5,7]))
    self.assertSetEqual(system.bonds[7],set([6]))
    self.assertSetEqual(system.bonds[8],set([9]))
    self.assertSetEqual(system.bonds[9],set([8,10]))
    self.assertSetEqual(system.bonds[10],set([9,11]))
    self.assertSetEqual(system.bonds[11],set([10]))

    system.remove_molecule(molecule=chain2,remove_beads=True)
    self.assertEqual(system.nbeads,len(x1+x3))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1],chain3)
    self.assertIs(system.molecules[1].system,system)
    self.assertSetEqual(system.molecules[0].indices,set(range(4)))
    self.assertSetEqual(system.molecules[1].indices,set(range(4,8)))
    self.assertListEqual(list(system.x),x1+x3)
    self.assertListEqual(list(system.y),y1+y3)
    self.assertListEqual(list(system.z),z1+z3)
    self.assertListEqual(list(system.types),types1+types3)
    self.assertSetEqual(system.bonds[0],set([1]))
    self.assertSetEqual(system.bonds[1],set([0,2]))
    self.assertSetEqual(system.bonds[2],set([1,3]))
    self.assertSetEqual(system.bonds[3],set([2]))
    self.assertSetEqual(system.bonds[4],set([5]))
    self.assertSetEqual(system.bonds[5],set([4,6]))
    self.assertSetEqual(system.bonds[6],set([5,7]))
    self.assertSetEqual(system.bonds[7],set([6]))



if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(Molecule_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


