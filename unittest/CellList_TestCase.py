import unittest
import numpy as np
from typySim.cy import CellList


class CellList_TestCase(unittest.TestCase):
  def create_cell_list(self,cell_grid,box):
    self.cellList = CellList(*cell_grid)
    self.cellList.set_box_size(*box)
  def create_dummy_positions(self,cell_grid,box,central_origin):
    grid_size = [b/float(d) for b,d in zip(box,cell_grid)]
    dx = cell_grid[0]
    dy = cell_grid[1]
    dz = cell_grid[2]
    bx = box[0]
    by = box[1]
    bz = box[2]
    gx = bx/float(dx)
    gy = by/float(dy)
    gz = bz/float(dz)
    position_list = []
    for ix in range(dx):
      for iy in range(dy):
        for iz in range(dz):
          position_list.append([ix*gx+gx/2.0,iy*gy+gy/2.0,iz*gz+gz/2.0])
    self.position_array = np.array(position_list)
    if central_origin:
      self.position_array[:,0] -= (bx/2.0)
      self.position_array[:,1] -= (by/2.0)
      self.position_array[:,2] -= (bz/2.0)
    return self.position_array.shape[0]
  def generate_random_position(self,box,central_origin):
    x = box[0]*np.random.random(1)
    y = box[1]*np.random.random(1)
    z = box[2]*np.random.random(1)
    if central_origin:
      x-=box[0]/2.0
      y-=box[1]/2.0
      z-=box[2]/2.0
    return x,y,z
  def base_box_test(self,box,cell_grid,central_origin):
    ncells = cell_grid[0]*cell_grid[1]*cell_grid[2]

    # Set up the cell list and positions array for this box
    self.create_cell_list(cell_grid=cell_grid,box=box)
    nbeads = self.create_dummy_positions(cell_grid=cell_grid,box=box,central_origin=central_origin)
    self.cellList.build_nlist(self.position_array,central_origin)

    #There should be 27 neighbors for each bead based on the dummy position array
    neighs = self.cellList.get_neighbors_by_pos(self.position_array[0])
    self.assertEqual(len(neighs),27)
    neighs = self.cellList.get_neighbors_by_pos(self.position_array[10])
    self.assertEqual(len(neighs),27)
    neighs = self.cellList.get_neighbors_by_pos(self.position_array[-1])
    self.assertEqual(len(neighs),27)

    # generate 100 random positions and add them to the CellList
    for i in range(100):
      x,y,z = self.generate_random_position(box,central_origin)
      self.cellList.insert_bead(nbeads+i,x,y,z)

    # remove all 100 random positions
    for i in range(100):
      self.cellList.remove_bead(nbeads+i)

    # test that the base CellList parameters are as expected or unchanged
    nlist = self.cellList.get_nlist()
    self.assertEqual(nlist['nbeads'],nbeads+100)
    self.assertEqual(nlist['ncells'],ncells)
    for i,key in enumerate(['bx','by','bz']):
      self.assertEqual(nlist[key],box[i])
    for i,key in enumerate(['nx','ny','nz']):
      self.assertEqual(nlist[key],cell_grid[i])
    for i,key in enumerate(['dx','dy','dz']):
      self.assertEqual(nlist[key],box[i]/float(cell_grid[i]))

    #There should be 27 neighbors for each bead based on the dummy position array
    neighs = self.cellList.get_neighbors_by_pos(self.position_array[12])
    self.assertEqual(len(neighs),27)

    # Does nlist resetting behave as expected?
    self.cellList.build_nlist(self.position_array,central_origin)
    nlist = self.cellList.get_nlist()
    self.assertEqual(nlist['nbeads'],nbeads)
    self.assertEqual(nlist['ncells'],ncells)
    for i,key in enumerate(['bx','by','bz']):
      self.assertEqual(nlist[key],box[i])
    for i,key in enumerate(['nx','ny','nz']):
      self.assertEqual(nlist[key],cell_grid[i])
    for i,key in enumerate(['dx','dy','dz']):
      self.assertEqual(nlist[key],box[i]/float(cell_grid[i]))

  def test_box1_N10x10x10_L50x50x50_central_origin(self):
    # Define the box parameters
    box = [50,50,50]
    cell_grid = [10,10,10]
    central_origin = True
    self.base_box_test(box,cell_grid,central_origin)


if __name__=="__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(CellList_TestCase)
  unittest.TextTestRunner(verbosity=2).run(suite)


