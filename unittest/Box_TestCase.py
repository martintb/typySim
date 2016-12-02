import unittest
import numpy as np
from typySim.core.cy.Box import Box


class Box_TestCase(unittest.TestCase):
  def test_box(self):

    box_lengths = [22,33,44.2]
    box = Box()

    # Default values
    self.assertTupleEqual(box.L,(-1,-1,-1))

    # L setter
    box.L = box_lengths[0]
    self.assertTupleEqual(box.L,(box_lengths[0],box_lengths[0],box_lengths[0]))

    # lx setter
    box.lx = -1
    self.assertTupleEqual(box.L,(-1,box_lengths[0],box_lengths[0]))

    # ly setter
    box.ly = box_lengths[1]
    self.assertTupleEqual(box.L,(-1,box_lengths[1],box_lengths[0]))

    # lz setter
    box.lz = box_lengths[2]
    # self.assertTupleEqual(box.L,(-1,box_lengths[1],box_lengths[2]))
    self.assertAlmostEqual(box.L[0],-1,delta=0.01)
    self.assertAlmostEqual(box.L[1],box_lengths[1],delta=0.01)
    self.assertAlmostEqual(box.L[2],box_lengths[2],delta=0.01)

    # half_L 
    # self.assertTupleEqual(box.half_L,(-1/2.0,box_lengths[1]/2.0,box_lengths[2]/2.0))
    self.assertAlmostEqual(box.half_L[0],-1/2.0,delta=0.01)
    self.assertAlmostEqual(box.half_L[1],box_lengths[1]/2.0,delta=0.01)
    self.assertAlmostEqual(box.half_L[2],box_lengths[2]/2.0,delta=0.01)

    #lo/hi
    self.assertAlmostEqual(box.xhi,-1/2.0,delta=0.01)
    self.assertAlmostEqual(box.xlo,1/2.0,delta=0.01)
    self.assertAlmostEqual(box.yhi,box_lengths[1]/2.0,delta=0.01)
    self.assertAlmostEqual(box.ylo,-box_lengths[1]/2.0,delta=0.01)
    self.assertAlmostEqual(box.zhi,box_lengths[2]/2.0,delta=0.01)
    self.assertAlmostEqual(box.zlo,-box_lengths[2]/2.0,delta=0.01)






if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(Box_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


