import random
import logging
import cPickle
import numpy as np

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
    self.name='MonteCarlo'
    self.rates = {}
    self.rates['total_attempted'] = 0
    self.rates['total_accepted'] = 0
    self.logger = logging.getLogger(__name__)
  def add_move(self,move):
    move.engine = self
    move.system = self.system
    self.moveList.append(move)
    if not (move.name in self.rates):
      self.rates[move.name] = 0
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
      viz.draw_all(bonds=False)
      viz.draw_box()
      viz.show(blocking=False)

    if pkl_rate is not None:
      pkl = {}



    self.TPE = self.system.get_compute('TotalPotentialEnergy')
    self.TPE_list = []
    self.TPE_list.append(sum(self.TPE.compute()))

    for i in range(num_attempts):
      move = random.choice(self.moveList)
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
        logStr +='Step {}/{}, rate: {:4.3f} U: {}'.format(i,num_attempts-1,rate,mc_move_data['string']) + bcolors.ENDC
        logStr += bcolors.ENDC
        self.logger.info(logStr)

      if (viz is not None) and success and (i%log_rate)==0:
        viz.clear()
        viz.draw_all(bonds=False)
        viz.draw_box()
        viz.show(blocking=False)

      if (pkl_rate is not None) and (i%pkl_rate)==0:
        pkl[i] = {}
        pkl[i]['x'] = np.array(self.system.x)
        pkl[i]['y'] = np.array(self.system.y)
        pkl[i]['z'] = np.array(self.system.z)
        pkl[i]['t'] = np.array(self.system.types)
        pkl[i]['L'] = np.array(self.system.box.L)
        pkl[i]['bonds'] = np.array(self.system.bonds.get_pairlist())
        pkl[i]['move_data'] = mc_move_data

    accepted = self.rates['total_accepted']
    attempted = self.rates['total_attempted']
    self.logger.info('Acceptance rate: {}/{} = {}'.format(accepted,attempted,rate))

    if pkl_rate is not None:
      self.logger.info('Logging trajectory to {}'.format(pkl_name))
      with open(pkl_name,'wb') as f:
        cPickle.dump(pkl,f,-1)
