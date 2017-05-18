import subprocess
import os
import copy
import logging
import numpy as np
import time
from .. import util

class Packmol(object):
  '''
  :class:`Engine` for conducting generating initial conditions for futher simulations. This
  Engine uses the Packmol library

  XXX DETAILS TO FOLLOW...

  '''
  def __init__(self):
    super(Packmol,self).__init__() #must call parent class' constructor
    self.name='MonteCarlo'
    self.mol = {}
    self.logger = logging.getLogger(__name__)
  def add_molecule(self,name,mol,molData,num):
    self.mol[name] = {'obj':mol,'size':len(molData['x']),'data':molData,'name':name,'num':num}
  def write_xyz(self,fname,mol):
    nbeads = mol['size']
    t = mol['data']['types']
    x = mol['data']['x']
    y = mol['data']['y']
    z = mol['data']['z']
    with open(fname,'w') as f:
      f.write('{}\n'.format(nbeads))
      f.write('typySim generated XYZ file for molecule: {}\n'.format(mol['name']))
      for i in range(nbeads):
        f.write('{:>3s} {:10.5f} {:10.5f} {:10.5f}\n'.format(str(t[i]),x[i],y[i],z[i]))
  def run(self,tol=2.0,seed=np.random.randint(999999),boxScale=1.0,randomize=True,discale=1.1):
    start_time = time.time()

    packStr = ''
    packStr += 'tolerance {}\n'.format(tol)
    packStr += 'filetype xyz\n'
    # packStr += 'add_box_sides 1.0\n'
    packStr += 'discale {} \n'.format(discale)
    packStr += 'output packed.xyz\n'
    packStr += 'seed {:d}\n'.format(seed)
    if randomize:
      packStr += 'randominitialpoint\n'
    packStr += '\n\n'

    boxDict = {}
    boxDict['xlo'] = self.system.box.xlo*boxScale
    boxDict['ylo'] = self.system.box.ylo*boxScale
    boxDict['zlo'] = self.system.box.zlo*boxScale
    boxDict['xhi'] = self.system.box.xhi*boxScale
    boxDict['yhi'] = self.system.box.yhi*boxScale
    boxDict['zhi'] = self.system.box.zhi*boxScale
    boxStr = 'inside box {xlo} {ylo} {zlo} {xhi} {yhi} {zhi}\n'.format(**boxDict)

    self.logger.info('Converting molecules to xyz files...')
    for i,(name,mol) in enumerate(self.mol.items(),start=1):
      fname = name + '.xyz'
      # self.write_xyz(fname,mol)
      util.write_xyz(fname,x=mol['data']['x'],y=mol['data']['y'],z=mol['data']['z'],types=mol['data']['types'])

      packStr += 'structure {}\n'.format(fname)
      packStr += '\tnumber {:d}\n'.format(int(mol['num']))
      packStr += '\t' + boxStr
      if 'radii' in mol['data']:
        for radius,idex_list in mol['data']['radii'].iteritems():
          idex_str = ' '.join(str(idex) for idex in idex_list)
          packStr += '\tatoms {}\n'.format(idex_str)
          packStr += '\t\tradius {}\n'.format(radius)
          packStr += '\tend atoms\n'
      packStr += 'end structure\n'.format(fname)
      packStr += '\n\n'


    self.logger.info('Writing PACKMOL input...')
    with open('packmol.in','w') as f:
      f.write(packStr)

    self.logger.info('Calling PACKMOL...')
    try:
      with open('packmol.in','r') as f:
        subprocess.check_call('packmol',stdin=f)
    except OSError as e:
      self.logger.error('Packmol failed. This is likely because Packmol isn\'t properly installed.')
      raise e


    if not os.path.isfile('packed.xyz'):
      self.logger.error('Packmol failed. See Packmol output for details.')
      raise RuntimeError()

    self.logger.info('Unpacking PACKMOL results...')
    packedPos = list(np.loadtxt('packed.xyz',skiprows=2))
    numMol = len(self.mol)
    for i,(name,mol) in enumerate(self.mol.items(),start=1):
      for j in range(mol['num']):
        molPos = []
        molTypes = []
        for n in range(mol['size']):
          t,x,y,z = packedPos.pop(0)
          molPos.append([x,y,z])
          molTypes.append(t)
        molPos = np.array(molPos)

        if not np.all(mol['data']['types'] == molTypes):
          raise ValueError('Non-matching types in system and packed.xyz!')

        molData = {}
        molData.update(mol['data'])
        molData['x'] = molPos[:,0]
        molData['y'] = molPos[:,1]
        molData['z'] = molPos[:,2]
        molData['bond_shift'] = True
        molObj = copy.deepcopy(mol['obj'])
        self.system.add_molecule(molObj,**molData)

    end_time = time.time()
    run_time_s = end_time - start_time
    run_time_h = (end_time - start_time)/3600
    self.logger.info('Packmol run finished!')
    self.logger.info('Total Execution Time: {} seconds or {} hours'.format(run_time_s,run_time_h))

