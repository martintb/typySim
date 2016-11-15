import random
import logging

class MonteCarlo(object):
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
    print 'Acceptance rate: {}/{} = {}'.format(accepted,attempted,rate)
