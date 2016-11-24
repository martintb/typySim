
if __name__ == '__main__':
  import unittest
  from CellList_TestCase import CellList_TestCase
  from Box_TestCase import Box_TestCase
  from System_TestCase import System_TestCase
  from Molecule_TestCase import Molecule_TestCase

  suite_list = []
  suite_list.append(unittest.TestLoader().loadTestsFromTestCase(CellList_TestCase))
  suite_list.append(unittest.TestLoader().loadTestsFromTestCase(Box_TestCase))
  suite_list.append(unittest.TestLoader().loadTestsFromTestCase(System_TestCase))
  suite_list.append(unittest.TestLoader().loadTestsFromTestCase(Molecule_TestCase))
  suite = unittest.TestSuite(suite_list)
  unittest.TextTestRunner(verbosity=2).run(suite)
