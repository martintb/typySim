from numpy.random import choice
import logging
import cPickle
import numpy as np

import ipdb; ist = ipdb.set_trace

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class MonteCarlo(object):
  '''
  Base :class:`Engine` for conducting  Monte Carlo (MC) style simulations

  This engine accumulates a list of MC moves which subclass :class:`MonteCarloMove`
  and then randomly chooses from them at every timestep. After calling the moves
  :func:`attempt` method the move records whether the move was accepted or rejected.

  '''
  def __init__(self):
    super(MonteCarlo,self).__init__() #must call parent class' constructor
    self.moveList = list()
    self.moveWeights = list()
    self.name='MonteCarlo'
    self.rates = {}
    self.rates['total_attempted'] = 0
    self.rates['total_accepted'] = 0
    self.viz = None
    self.logger = logging.getLogger(__name__)
  def add_move(self,move,weight=1):
    move.set_engine(self)
    self.moveList.append(move)
    self.moveWeights.append(weight)
    if not (move.name in self.rates):
      self.rates[move.name] = 0
  def remove_move(self,move):
    self.moveList.remove(move)
    del self.rates[move.name]
    del move
  def run(self,num_attempts,log_rate = 5,viz=None,pkl_rate=None,pkl_name='trj.pkl'):
    '''
    Contduct the Monte Carlo simulation.

    Parameters
    ----------
    num_attempts : int, *Required*
        The number of MC moves that the engine should attempt in succession before
        stopping. 
    log_rate : int, *Optional*
        The rate at which the logger (at the INFO level) should log updates.
    '''
    # if self.system.neighbor_list is not None:
    # self.system.neighbor_list.build_nlist(self.system.x,self.s,central_origin=True)

    if viz is not None:
      self.viz = viz
      viz.draw_system(bonds=True)
      viz.show(blocking=False,resetCamera=False)

    if pkl_rate is not None:
      pkl = {}

    self.TPE = self.system.get_compute('TotalPotentialEnergy')
    self.TPE_list = []
    self.TPE_list.append(sum(self.TPE.compute()))

    totalWeight = float(sum(self.moveWeights))
    move_probs = [i/totalWeight for i in self.moveWeights]
    for i in range(num_attempts):
      move = choice(self.moveList,p=move_probs)
      success,mc_move_data = move.attempt()

      self.rates['total_attempted'] += 1
      if success:
        self.rates['total_accepted'] += 1
        self.rates[move.name] += 1
        self.TPE_list.append(mc_move_data['U'])
      else:
        self.TPE_list.append(self.TPE_list[-1])
      if (i%log_rate)==0:
        accepted = self.rates['total_accepted']
        attempted = self.rates['total_attempted']
        rate = accepted/float(attempted)
        if success:
          color = bcolors.OKGREEN
        else:
          color = bcolors.FAIL
        logStr  = color 
        # logStr +='Step {}/{}, rate: {:4.3f} U: {} W: {}'.format(i,num_attempts-1,rate,self.TPE_list[-1],mc_move_data['string']) + bcolors.ENDC
        logStr +='Step {}/{}, rate: {:4.3f} {}'.format(i,num_attempts-1,rate,mc_move_data['string']) + bcolors.ENDC
        logStr += bcolors.ENDC
        self.logger.info(logStr)

      if (viz is not None) and success and (i%log_rate)==0:
        viz.clear()
        # viz.draw_system(bonds=False)
        viz.draw_system(bonds=True)
        viz.show(blocking=False,resetCamera=False)

      if (pkl_rate is not None) and (i%pkl_rate)==0:
        pkl[i] = {}
        pkl[i]['nbeads']    = np.array(self.system.nbeads)
        pkl[i]['L']         = np.array(self.system.box.L)
        pkl[i]['move_data'] = mc_move_data
        pkl[i]['x']         = np.array(self.system.x)
        pkl[i]['y']         = np.array(self.system.y)
        pkl[i]['z']         = np.array(self.system.z)
        pkl[i]['imx']       = np.array(self.system.imx)
        pkl[i]['imy']       = np.array(self.system.imy)
        pkl[i]['imz']       = np.array(self.system.imz)
        pkl[i]['types']     = np.array(self.system.types)
        pkl[i]['bonds']     = np.array(self.system.bonds.get_pairlist())
        pkl[i]['molecules'] = []
        for mol in self.system.molecules:
          mol_dict = {}
          mol_dict['name']       = mol.name
          mol_dict['indices']    = mol.indices
          mol_dict['properties'] = mol.properties
          pkl[i]['molecules'].append(mol_dict)


    accepted = self.rates['total_accepted']
    attempted = self.rates['total_attempted']
    self.logger.info('Acceptance rate: {}/{} = {}'.format(accepted,attempted,rate))

    if pkl_rate is not None:
      self.logger.log(logging.INFO+1,'Logging trajectory to {}'.format(pkl_name))
      with open(pkl_name,'wb') as f:
        cPickle.dump(pkl,f,-1)

    if (viz is not None):
      viz.clear()
      viz.draw_system(bonds=True)
      viz.show()
