import cPickle
import numpy as np
from numpy.random import randint,random,choice

from typySim import molecule
from typySim.geometry import linalg
from typySim.engine.MonteCarloMove import MonteCarloMove
from typySim.engine.RosenbluthChain import RosenbluthChain

import ipdb as pdb; ist = pdb.set_trace


class CBMCMove(MonteCarloMove):
  '''
  Definitions
  -----------
  local index :  [0,len(mol.indices))
    Integer ranging between 0 and the length of the chain being regrown. Used
    to access the indices and data arrays for the regrowing chain.

  sys index : [0,self.system.nbeads)
    Global index of a bead in the system arrays (e.g. self.system.x ). This 
    value ranges from .

  new index : [self.system.nbeads,self.system.nbeads+len(mol.indices))
    Global index of the trial_move bead. This value ranges from self.system to
    the number of beads being regrown

  XXX_index_trial : [0,self.num_trials)
    Represents the index of one of the possibi

  '''
  def __init__(self,regrowth_types,num_trials=25,regrowth_min=2,regrowth_max=1000):
    super(CBMCMove,self).__init__() #must call parent class' constructor
    self.name='CBMCMove'
    self.regrowth_types = regrowth_types

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max

    self.bias = None
    self.rosen_chain = RosenbluthChain(num_trials,regrowth_min,regrowth_max,self.bias)
  def set_engine(self,engine):
    self.engine               = engine
    self.system               = engine.system
    self.rosen_chain.engine        = engine
    self.rosen_chain.system        = engine.system
  @MonteCarloMove.counter
  def attempt(self):
    mc_move_data = {}

    ################################
    ## DETERMINE REGROWTH INDICES##
    ################################
    # Choose a chain molecule and a starting point
    mol = self.system.select.random_molecule(names=['ChainSegment'])
    indices = list(mol.indices)
    chain_length = len(indices)

    if chain_length<self.regrowth_min+2:
      accept = False
      mc_move_data['string'] = 'chain_is_too_short'
      return accept,mc_move_data

    self.rosen_chain.set_indices(mol.indices,internal_growth=False)
    self.rosen_chain.build_arrays()
    self.rosen_chain.copy_retrace_to_regrowth()
    self.rosen_chain.acceptance_package['bonds'] = {'added':[], 'removed':[]}
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[],'added':[],'removed':[]}

    ################################
    ## INITIAL ENERGY CALCULATION ##
    ################################
    # Need to calculate the base energy of the system that we are regrowing into
    UOld    = self.engine.TPE_list[-1]
    URegrow = self.engine.TPE.compute(partial_indices=self.rosen_chain.indices)
    UBase   = UOld-sum(URegrow)

    ####################
    ## GROW NEW CHAIN ##
    ####################
    abort,rosen_weights_new,bias_weights_new,Jnew = self.rosen_chain.calc_rosenbluth(UBase,retrace=False)
    if abort:
      mc_move_data['string'] = 'bad_trial_move_new'
      accept = False
      return accept,mc_move_data

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    abort,rosen_weights_old,bias_weights_old,Jold = self.rosen_chain.calc_rosenbluth(UBase,retrace=True)
    if abort:
      mc_move_data['string'] = 'bad_trial_move_old'
      accept = False
      return accept,mc_move_data

    #################################
    ## CALCULATE FULL ROSEN FACTOR ##
    #################################
    # print 'TOPROSEN',np.sum(rosen_weights_new,axis=1)
    # print 'TOPBIAS',bias_weights_old
    # print 'BOTBIAS',np.sum(rosen_weights_old,axis=1)
    # print 'BOTROSEN',bias_weights_new
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))*np.product(bias_weights_old)*Jnew
    Wold = np.product(np.sum(rosen_weights_old,axis=1))*np.product(bias_weights_new)*Jold
    mc_move_data['Wnew'] = Wnew
    mc_move_data['Wold'] = Wold
    mc_move_data['k'] = self.rosen_chain.length

    #######################
    ## ACCEPT OR REJECT? ##
    #######################
    if Wnew>Wold:
      accept = True
    else:
      ranf = np.random.rand()
      if ranf<(Wnew/Wold):
        accept=True
      else:
        accept=False

    if accept:
      mc_move_data['U'] = self.rosen_chain.acceptance_package['U']
      self.rosen_chain.apply_acceptance_package()

    mc_move_data['string'] = 'CBMC::Wnew/Wold: {:3.2e}/{:3.2e}={:5.4f}'.format(Wnew,Wold,Wnew/Wold)
    return accept,mc_move_data









