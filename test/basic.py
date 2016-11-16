import numpy as np
import typySim
import pdb; st = pdb.set_trace
import logging,sys
logging.basicConfig(stream=sys.stdout,level=logging.DEBUG)


sys = typySim.core.System(cell_grid=(5,2,10))
sys.box.L=10
sys.box.lx=25
sys.box.lz=50

hexp_kwargs = {}
hexp_kwargs['lx'] = sys.box.lx
hexp_kwargs['ly'] = sys.box.ly
hexp_kwargs['nz'] = 5
hexp_kwargs['diameter'] = 1.0
hexp_kwargs['topType'] = 1
hexp_kwargs['middleType'] = 0
hexp_kwargs['bottomType'] = 1
hexp1 = typySim.molecule.HexagonalSurface()
molData,boxData = hexp1.build(**hexp_kwargs)
new_types = []
for t in molData['types']:
  if t==0:
    new_types.append(t)
  elif (t!=0) and (np.random.random()>0.5):
    new_types.append(1)
  else:
    new_types.append(0)
molData['types'] = new_types
sys.add_molecule(hexp1,molData = molData)
sys.lx = boxData['lx']
sys.ly = boxData['ly']


length = 8
ChainSeg = typySim.molecule.ChainSegment()
minZ = 2
maxZ = 8
dz = float(maxZ-minZ)/float(length)
pos = 10*np.random.random((length,3))-5
pos[:,2] = np.arange(minZ,maxZ,dz)-5
molData['positions'] = pos
molData['types'] = [4]*length
molData['bonds'] = [[i,j] for i,j in zip(range(length-1),range(1,length))]
# sys.add_molecule(ChainSeg,molData=molData)
sys.reset_all_molecules()

sys.molecules[0].positions[:,2]+=((sys.box.lz/2.0) - 0.5)
sys.box.wrap_all_positions()

sys.add_engine(typySim.engine.MonteCarlo())
RWGM = typySim.engine.RandomWalkGrowthMove()
RWGM.grow_from_types = [1,4]
RWGM.chain_end_type = 4
RWGM.chain_middle_type = 3
sys.engine.add_move(RWGM)
sys.engine.run(100000,log_rate=100)

viz1 = typySim.viz.SystemViewer(sys)
# viz1.add_system(check_bonds=True)
viz1.draw_all()
viz1.draw_box()
viz1.show()




