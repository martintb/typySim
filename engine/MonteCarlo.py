import random
import logging

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
  def run(self,num_attempts,log_rate = 5):
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
    if self.system.box.cellList is not None:
      self.system.box.cellList.build_nlist(self.system.positions,central_origin=True)

    for i in range(num_attempts):
      move = random.choice(self.moveList)
      success = move.attempt()

      self.rates['total_attempted'] += 1
      if success:
        self.rates['total_accepted'] += 1
        self.rates[move.name] += 1
      if (i%log_rate)==0:
        accepted = self.rates['total_accepted']
        attempted = self.rates['total_attempted']
        rate = accepted/float(attempted)
        logStr  = 'Step {}/{}, rate: {}'.format(i,num_attempts-1,rate)
        self.logger.info(logStr)

    accepted = self.rates['total_accepted']
    attempted = self.rates['total_attempted']
    rate = accepted/float(attempted)
    self.logger.info(print 'Acceptance rate: {}/{} = {}'.format(accepted,attempted,rate))
