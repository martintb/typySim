import unittest
import numpy as np
from typySim.core import Box


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
    self.assertTupleEqual(box.L,(-1,box_lengths[1],box_lengths[2]))

    # half_L 
    self.assertTupleEqual(box.half_L,(-1/2.0,box_lengths[1]/2.0,box_lengths[2]/2.0))

    #lo/hi
    t1 = (box.xlo,box.ylo,box.zlo)
    t2 = (box.xhi,box.yhi,box.zhi)
    t3 = (1/2.0,-box_lengths[1]/2.0,-box_lengths[2]/2.0)
    t4 = (-1/2.0,box_lengths[1]/2.0,box_lengths[2]/2.0)
    self.assertTupleEqual(t1,t3)
    self.assertTupleEqual(t2,t4)






if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(Box_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


