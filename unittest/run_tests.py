
if __name__ == '__main__':
  import unittest
  from CellList_TestCase import CellList_TestCase
  from Box_TestCase import Box_TestCase

  suite_list = []
  suite_list.append(unittest.TestLoader().loadTestsFromTestCase(CellList_TestCase))
  suite_list.append(unittest.TestLoader().loadTestsFromTestCase(Box_TestCase))
  suite = unittest.TestSuite(suite_list)
  unittest.TextTestRunner(verbosity=2).run(suite)
