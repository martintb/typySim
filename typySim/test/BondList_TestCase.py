import unittest
import numpy as np
from typySim.core import BondList


class BondList_TestCase(unittest.TestCase):
  def test_add_remove_bonds(self):
    BL = BondList()

    BL.expand(10)
    BL.add(4,5,0)
    BL.add(4,6,0)
    BL.add(4,0,0)

    np.testing.assert_array_equal(BL[0],np.array([ 4,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[1],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[2],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[3],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[4],np.array([ 0, 5, 6,-1,-1]))
    np.testing.assert_array_equal(BL[5],np.array([ 4,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[6],np.array([ 4,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[7],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[8],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[9],np.array([-1,-1,-1,-1,-1]))

    BL.remove(4,5,0)

    np.testing.assert_array_equal(BL[0],np.array([ 4,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[1],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[2],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[3],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[4],np.array([ 0, 6,-1,-1,-1]))
    np.testing.assert_array_equal(BL[5],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[6],np.array([ 4,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[7],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[8],np.array([-1,-1,-1,-1,-1]))
    np.testing.assert_array_equal(BL[9],np.array([-1,-1,-1,-1,-1]))




if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(BondList_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


