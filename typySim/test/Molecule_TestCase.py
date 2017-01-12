import unittest
import numpy as np
from typySim.core import System
from typySim import molecule
import ipdb; ist = ipdb.set_trace


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
    system.add_molecule(chain,x=x,y=y,z=z,types=types,bonds=bonds,bond_shift=True)


    self.assertEqual(system.nbeads,len(x))
    self.assertIs(system.molecules[0],chain)
    self.assertIs(system.molecules[0].system,system)
    self.assertListEqual(system.molecules[0].indices,range(8))
    self.assertListEqual(list(system.x),x)
    self.assertListEqual(list(system.y),y)
    self.assertListEqual(list(system.z),z)
    self.assertListEqual(list(chain.x.compressed()),x)
    self.assertListEqual(list(chain.y.compressed()),y)
    self.assertListEqual(list(chain.z.compressed()),z)
    self.assertListEqual(list(chain.types.compressed()),types)
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2, 4,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[4],np.array([3, 5,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[5],np.array([4, 6,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[6],np.array([5, 7,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[7],np.array([6,-1,-1,-1,-1]))
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
    system.add_molecule(chain1,x=x1,y=y1,z=z1,types=types1,bonds=bonds1,bond_shift=True)
    system.add_molecule(chain2,x=x2,y=y2,z=z2,types=types2,bonds=bonds2,bond_shift=True)


    self.assertEqual(system.nbeads,len(x1+x2))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[1],chain2)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1].system,system)
    self.assertListEqual(system.molecules[0].indices,range(4))
    self.assertListEqual(system.molecules[1].indices,range(4,8))
    self.assertListEqual(list(system.x),x1+x2)
    self.assertListEqual(list(system.y),y1+y2)
    self.assertListEqual(list(system.z),z1+z2)
    self.assertListEqual(list(system.types),types1+types2)
    self.assertListEqual(list(chain1.x.compressed()),x1)
    self.assertListEqual(list(chain1.y.compressed()),y1)
    self.assertListEqual(list(chain1.z.compressed()),z1)
    self.assertListEqual(list(chain1.types.compressed()),types1)
    self.assertListEqual(list(chain2.x.compressed()),x2)
    self.assertListEqual(list(chain2.y.compressed()),y2)
    self.assertListEqual(list(chain2.z.compressed()),z2)
    self.assertListEqual(list(chain2.types.compressed()),types2)
    # self.assertSetEqual(system.bonds[0],set([1]))
    # self.assertSetEqual(system.bonds[1],set([0,2,4]))
    # self.assertSetEqual(system.bonds[2],set([1,3]))
    # self.assertSetEqual(system.bonds[3],set([2]))
    # self.assertSetEqual(system.bonds[4],set([1,5]))
    # self.assertSetEqual(system.bonds[5],set([4,6]))
    # self.assertSetEqual(system.bonds[6],set([5,7]))
    # self.assertSetEqual(system.bonds[7],set([6]))
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2, 4,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[4],np.array([1, 5,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[5],np.array([4, 6,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[6],np.array([5, 7,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[7],np.array([6,-1,-1,-1,-1]))

    system.remove_molecule(molecule=chain1,remove_beads=True)
    self.assertEqual(system.nbeads,len(x2))
    self.assertIs(system.molecules[0],chain2)
    self.assertIs(system.molecules[0].system,system)
    self.assertListEqual(system.molecules[0].indices,range(4))
    self.assertListEqual(list(system.x),x2)
    self.assertListEqual(list(system.y),y2)
    self.assertListEqual(list(system.z),z2)
    self.assertListEqual(list(system.types),types2)
    self.assertListEqual(list(chain2.x.compressed()),x2)
    self.assertListEqual(list(chain2.y.compressed()),y2)
    self.assertListEqual(list(chain2.z.compressed()),z2)
    self.assertListEqual(list(chain2.types.compressed()),types2)
    # self.assertSetEqual(system.bonds[0],set([1]))
    # self.assertSetEqual(system.bonds[1],set([0,2]))
    # self.assertSetEqual(system.bonds[2],set([1,3]))
    # self.assertSetEqual(system.bonds[3],set([2]))
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2,-1,-1,-1,-1]))
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
    system.add_molecule(chain1,x=x1,y=y1,z=z1,types=types1,bonds=bonds1,bond_shift=True)
    system.add_molecule(chain2,x=x2,y=y2,z=z2,types=types2,bonds=bonds2,bond_shift=True)
    system.add_molecule(chain3,x=x3,y=y3,z=z3,types=types3,bonds=bonds3,bond_shift=True)



    self.assertEqual(system.nbeads,len(x1+x2+x3))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[1],chain2)
    self.assertIs(system.molecules[2],chain3)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1].system,system)
    self.assertIs(system.molecules[2].system,system)
    self.assertListEqual(system.molecules[0].indices,range(4))
    self.assertListEqual(system.molecules[1].indices,range(4,8))
    self.assertListEqual(system.molecules[2].indices,range(8,12))
    self.assertListEqual(list(system.x),x1+x2+x3)
    self.assertListEqual(list(system.y),y1+y2+y3)
    self.assertListEqual(list(system.z),z1+z2+z3)
    self.assertListEqual(list(system.types),types1+types2+types3)
    # self.assertSetEqual(system.bonds[0],set([1]))
    # self.assertSetEqual(system.bonds[1],set([0,2,4]))
    # self.assertSetEqual(system.bonds[2],set([1,3]))
    # self.assertSetEqual(system.bonds[3],set([2]))
    # self.assertSetEqual(system.bonds[4],set([1,5]))
    # self.assertSetEqual(system.bonds[5],set([4,6]))
    # self.assertSetEqual(system.bonds[6],set([5,7]))
    # self.assertSetEqual(system.bonds[7],set([6]))
    # self.assertSetEqual(system.bonds[8],set([9]))
    # self.assertSetEqual(system.bonds[9],set([8,10]))
    # self.assertSetEqual(system.bonds[10],set([9,11]))
    # self.assertSetEqual(system.bonds[11],set([10]))
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2, 4,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[4],np.array([1, 5,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[5],np.array([4, 6,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[6],np.array([5, 7,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[7],np.array([6,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[8],np.array([9,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[9],np.array([8,10,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[10],np.array([9,11,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[11],np.array([10,-1,-1,-1,-1]))

    system.remove_molecule(molecule=chain2,remove_beads=True)
    self.assertEqual(system.nbeads,len(x1+x3))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1],chain3)
    self.assertIs(system.molecules[1].system,system)
    self.assertListEqual(system.molecules[0].indices,range(4))
    self.assertListEqual(system.molecules[1].indices,range(4,8))
    self.assertListEqual(list(system.x),x1+x3)
    self.assertListEqual(list(system.y),y1+y3)
    self.assertListEqual(list(system.z),z1+z3)
    self.assertListEqual(list(system.types),types1+types3)
    # self.assertSetEqual(system.bonds[0],set([1]))
    # self.assertSetEqual(system.bonds[1],set([0,2]))
    # self.assertSetEqual(system.bonds[2],set([1,3]))
    # self.assertSetEqual(system.bonds[3],set([2]))
    # self.assertSetEqual(system.bonds[4],set([5]))
    # self.assertSetEqual(system.bonds[5],set([4,6]))
    # self.assertSetEqual(system.bonds[6],set([5,7]))
    # self.assertSetEqual(system.bonds[7],set([6]))
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[4],np.array([5,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[5],np.array([4, 6,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[6],np.array([5, 7,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[7],np.array([6,-1,-1,-1,-1]))

  def test_add_indices(self):
    ''' 
    Add three chains to a system, remove the middle chain, and then 
    add those center beads to the third chain.
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
    system.add_molecule(chain1,x=x1,y=y1,z=z1,types=types1,bonds=bonds1,bond_shift=True)
    system.add_molecule(chain2,x=x2,y=y2,z=z2,types=types2,bonds=bonds2,bond_shift=True)
    system.add_molecule(chain3,x=x3,y=y3,z=z3,types=types3,bonds=bonds3,bond_shift=True)

    #remove chain2
    chain2_removed = system.remove_molecule(molecule=chain2,remove_beads=False)

    #add beads from chain2 to chain3
    chain3.add_indices(chain2_removed.indices)

    self.assertIs(chain2,chain2_removed)
    self.assertEqual(chain3.size,8)
    self.assertEqual(system.molecules[1].size,8)

    self.assertEqual(system.nbeads,len(x1+x2+x3))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[1],chain3)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1].system,system)
    self.assertListEqual(system.molecules[0].indices,range(4))
    self.assertListEqual(system.molecules[1].indices,range(4,12))
    self.assertListEqual(list(system.x),x1+x2+x3)
    self.assertListEqual(list(system.y),y1+y2+y3)
    self.assertListEqual(list(system.z),z1+z2+z3)
    self.assertListEqual(list(system.types),types1+types2+types3)
    self.assertListEqual(list(chain1.x.compressed()),x1)
    self.assertListEqual(list(chain1.y.compressed()),y1)
    self.assertListEqual(list(chain1.z.compressed()),z1)
    self.assertListEqual(list(chain1.types.compressed()),types1)
    self.assertListEqual(list(chain3.x.compressed()),x2+x3) #combined molecule1 and molecule2
    self.assertListEqual(list(chain3.y.compressed()),y2+y3) #combined molecule1 and molecule2
    self.assertListEqual(list(chain3.z.compressed()),z2+z3) #combined molecule1 and molecule2
    self.assertListEqual(list(chain3.types.compressed()),types2+types3) #combined molecule1 and molecule2
    # self.assertSetEqual(system.bonds[0],set([1]))
    # self.assertSetEqual(system.bonds[1],set([0,2,4]))
    # self.assertSetEqual(system.bonds[2],set([1,3]))
    # self.assertSetEqual(system.bonds[3],set([2]))
    # self.assertSetEqual(system.bonds[4],set([1,5]))
    # self.assertSetEqual(system.bonds[5],set([4,6]))
    # self.assertSetEqual(system.bonds[6],set([5,7]))
    # self.assertSetEqual(system.bonds[7],set([6]))
    # self.assertSetEqual(system.bonds[8],set([9]))
    # self.assertSetEqual(system.bonds[9],set([8,10]))
    # self.assertSetEqual(system.bonds[10],set([9,11]))
    # self.assertSetEqual(system.bonds[11],set([10]))
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2, 4,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[4],np.array([1, 5,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[5],np.array([4, 6,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[6],np.array([5, 7,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[7],np.array([6,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[8],np.array([9,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[9],np.array([8,10,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[10],np.array([9,11,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[11],np.array([10,-1,-1,-1,-1]))

    for mol in system.molecule_map:
      self.assertIsNot(mol,system.DummyMolecule)
  def test_distribute(self):
    ''' 
    Add three chains to a system, remove the middle chain, and then 
    add those center beads to the third chain.
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
    x3 = range(20,24)
    y3 = range(24,28)
    z3 = range(28,32)
    types2= [2,3,3,2]
    bonds2 = [[i,j] for i,j in zip(range(len(x2)-1),range(1,len(x2)))]
    bonds2 += [['1',0]] 
    types3= [3,4,4,3]
    bonds3 = [[i+4,j+4] for i,j in zip(range(len(x3)-1),range(1,len(x3)))]

    x23 = x2 + x3
    y23 = y2 + y3
    z23 = z2 + z3
    types23 = types2 + types3
    bonds23 = bonds2 + bonds3
    chain23 = molecule.ChainSegment()

    system = System()
    system.add_molecule(chain1,x=x1,y=y1,z=z1,types=types1,bonds=bonds1,bond_shift=True)
    system.add_molecule(chain23,x=x23,y=y23,z=z23,types=types23,bonds=bonds23,bond_shift=True)

    chain2, chain3 =  chain23.distribute([range(4,8),range(8,12)])
    
    # The new molecules should be their own references
    self.assertIsNot(chain2,chain23)
    self.assertIsNot(chain3,chain23)

    # ..but they should be the same type
    self.assertIs(type(chain2),type(chain23))
    self.assertIs(type(chain3),type(chain23))

    self.assertEqual(system.nbeads,len(x1+x2+x3))
    self.assertIs(system.molecules[0],chain1)
    self.assertIs(system.molecules[1],chain2)
    self.assertIs(system.molecules[2],chain3)
    self.assertIs(system.molecules[0].system,system)
    self.assertIs(system.molecules[1].system,system)
    self.assertIs(system.molecules[2].system,system)
    self.assertListEqual(system.molecules[0].indices,range(4))
    self.assertListEqual(system.molecules[1].indices,range(4,8))
    self.assertListEqual(system.molecules[2].indices,range(8,12))
    self.assertListEqual(list(system.x),x1+x2+x3)
    self.assertListEqual(list(system.y),y1+y2+y3)
    self.assertListEqual(list(system.z),z1+z2+z3)
    self.assertListEqual(list(system.types),types1+types2+types3)
    self.assertListEqual(list(chain1.x.compressed()),x1)
    self.assertListEqual(list(chain1.y.compressed()),y1)
    self.assertListEqual(list(chain1.z.compressed()),z1)
    self.assertListEqual(list(chain1.types.compressed()),types1)
    self.assertListEqual(list(chain2.x.compressed()),x2)
    self.assertListEqual(list(chain2.y.compressed()),y2)
    self.assertListEqual(list(chain2.z.compressed()),z2)
    self.assertListEqual(list(chain2.types.compressed()),types2)
    self.assertListEqual(list(chain3.x.compressed()),x3)
    self.assertListEqual(list(chain3.y.compressed()),y3)
    self.assertListEqual(list(chain3.z.compressed()),z3)
    self.assertListEqual(list(chain3.types.compressed()),types3)
    # self.assertSetEqual(system.bonds[0],set([1]))
    # self.assertSetEqual(system.bonds[1],set([0,2,4]))
    # self.assertSetEqual(system.bonds[2],set([1,3]))
    # self.assertSetEqual(system.bonds[3],set([2]))
    # self.assertSetEqual(system.bonds[4],set([1,5]))
    # self.assertSetEqual(system.bonds[5],set([4,6]))
    # self.assertSetEqual(system.bonds[6],set([5,7]))
    # self.assertSetEqual(system.bonds[7],set([6]))
    # self.assertSetEqual(system.bonds[8],set([9]))
    # self.assertSetEqual(system.bonds[9],set([8,10]))
    # self.assertSetEqual(system.bonds[10],set([9,11]))
    # self.assertSetEqual(system.bonds[11],set([10]))
    np.testing.assert_array_equal(system.bonds[0],np.array([1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[1],np.array([0, 2, 4,-1,-1]))
    np.testing.assert_array_equal(system.bonds[2],np.array([1, 3,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[3],np.array([2,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[4],np.array([1, 5,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[5],np.array([4, 6,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[6],np.array([5, 7,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[7],np.array([6,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[8],np.array([9,-1,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[9],np.array([8,10,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[10],np.array([9,11,-1,-1,-1]))
    np.testing.assert_array_equal(system.bonds[11],np.array([10,-1,-1,-1,-1]))

    for mol in system.molecule_map:
      self.assertIsNot(mol,system.DummyMolecule)
      self.assertIsNot(mol,chain23)







if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(Molecule_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


