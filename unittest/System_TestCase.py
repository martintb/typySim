import unittest
import numpy as np
from typySim.core import System


class System_TestCase(unittest.TestCase):
  def test_add_beads(self):
    system = System()

    #Basic sanity check that beads are being added
    x = [1,2,3,4]
    y = [3,4,5,6]
    z = [7,8,8,10]
    types = [0,1,1,0]
    bonds = [[0,1],[2,3],[0,1],[3,0]]
    system.add_beads(x=x,y=y,z=z,types=types,bonds=bonds)
    self.assertEqual(system.nbeads,len(x))
    self.assertListEqual(list(system.x),x)
    self.assertListEqual(list(system.y),y)
    self.assertListEqual(list(system.z),z)
    self.assertListEqual(list(system.types),types)
    self.assertSetEqual(system.bonds[0],set([1,3]))
    self.assertSetEqual(system.bonds[3],set([0,2]))

    #Check to make sure that, non-matching bead arrays are rejected
    x2 = [1]
    y2 = [3,4]
    z2 = [0]
    types2 = [0,1,3]
    self.assertRaises(ValueError,system.add_beads,x2,y2,z2,types2)

    bonds = [['0','5']]
    system.add_beads(x=x,y=y,z=z,types=types,bonds=bonds)
    self.assertEqual(system.nbeads,len(x*2))
    self.assertListEqual(list(system.x),x*2)
    self.assertListEqual(list(system.y),y*2)
    self.assertListEqual(list(system.z),z*2)
    self.assertListEqual(list(system.types),types*2)
    self.assertSetEqual(system.bonds[0],set([1,3,5]))
    self.assertSetEqual(system.bonds[3],set([0,2]))
    self.assertSetEqual(system.bonds[5],set([0]))
  def test_remove_beads(self):
    '''
    The trickiest part of adding and removing beads is correctly
    mapping the bondlists. Here is a description of this test
    is testing for. 

    Starting Polymer:
    0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 

    Remove 0,4,7:
    x - 1 - 2 - 3 - x - 5 - 6 - x 

        |   |   |       |   |
        V   V   V       V   V

        1 - 2 - 3       5 - 6     

    Renumbered:
        0 - 1 - 2       3 - 4     

        BONDLISTS
      OLD       NEW
    -------   -------
    0 : 1     0 : 1
    1 : 0,2   1 : 0,2
    2 : 1,3   2 : 1  
    3 : 2,4   3 : 4
    4 : 3,5   4 : 3
    5 : 4,6          
    6 : 5,7          
    7 : 6            
    '''
    system = System()

    #Basic sanity check that beads are being added
    x = range(8)
    y = range(8,16)
    z = range(16,24)
    types = [0,1,1,0]*2
    bonds = [[i,j] for i,j in zip(range(len(x)-1),range(1,len(x)))]
    system.add_beads(x=x,y=y,z=z,types=types,bonds=bonds)

    to_be_removed = [0,4,7]
    system.remove_beads(to_be_removed)
    x     = list(np.delete(x,to_be_removed))
    y     = list(np.delete(y,to_be_removed))
    z     = list(np.delete(z,to_be_removed))
    types = list(np.delete(types,to_be_removed))

    self.assertEqual(system.nbeads,len(x))
    self.assertListEqual(list(system.x),x)
    self.assertListEqual(list(system.y),y)
    self.assertListEqual(list(system.z),z)
    self.assertListEqual(list(system.types),types)
    self.assertSetEqual(system.bonds[0],set([1]))
    self.assertSetEqual(system.bonds[1],set([0,2]))
    self.assertSetEqual(system.bonds[2],set([1]))
    self.assertSetEqual(system.bonds[3],set([4]))
    self.assertSetEqual(system.bonds[4],set([3]))



if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(System_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


