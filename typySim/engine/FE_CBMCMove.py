import cPickle
import numpy as np
from numpy.random import randint,random,choice

from typySim import molecule
from typySim.geometry import linalg
from typySim.engine.MonteCarloMove import MonteCarloMove
from typySim.engine.RosenbluthChain import RosenbluthChain

import ipdb as pdb; ist = pdb.set_trace


class FE_CBMCMove(MonteCarloMove):
  '''
  Definitions
  -----------
  local index :  [0,len(regrowth_mol.indices))
    Integer ranging between 0 and the length of the chain being regrown. Used
    to access the indices and data arrays for the regrowing chain.

  sys index : [0,self.system.nbeads)
    Global index of a bead in the system arrays (e.g. self.system.x ). This 
    value ranges from .

  new index : [self.system.nbeads,self.system.nbeads+len(regrowth_mol.indices))
    Global index of the trial_move bead. This value ranges from self.system to
    the number of beads being regrown

  XXX_index_trial : [0,self.num_trials)
    Represents the index of one of the possibi

  '''
  def __init__(self,regrowth_types,bias_pkl,num_trials=25,regrowth_min=2,regrowth_max=6):
    super(FE_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='FE_CBMCMove'
    self.regrowth_types = regrowth_types

    with open(bias_pkl,'rb') as f:
      self.bias = cPickle.load(f)

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max
    self.chain1 = RosenbluthChain(num_trials,regrowth_min,regrowth_max,self.bias)
  def set_engine(self,engine):
    self.engine               = engine
    self.system               = engine.system
    self.chain1.engine        = engine
    self.chain1.system        = engine.system
  @MonteCarloMove.counter
  def attempt(self):
    mc_move_data = {}

    ################################
    ## DETERMINE REGROWTH INDICES##
    ################################
    # Choose a chain molecule and a starting point
    regrowth_mol = self.system.select.random_molecule(names=['ChainSegment'])
    indices = list(regrowth_mol.indices)
    chain_length = len(indices)

    if chain_length<self.regrowth_min+2:
      accept = False
      mc_move_data['string'] = 'chain_is_too_short'
      return accept,mc_move_data

    self.chain1.set_indices(regrowth_mol,internal_growth=True)
    self.chain1.build_old_arrays()
    self.chain1.new_bonds   = self.chain1.old_bonds
    self.chain1.new_anchors = self.chain1.old_anchors
    self.chain1.new_beacons = self.chain1.old_beacons

    ################################
    ## INITIAL ENERGY CALCULATION ##
    ################################
    # Need to calculate the base energy of the system that we are regrowing into
    UOld    = self.engine.TPE_list[-1]
    URegrow = self.engine.TPE.compute(partial_indices=self.chain1.regrowth_indices)
    UBase   = UOld-sum(URegrow)

    ####################
    ## GROW NEW CHAIN ##
    ####################
    rosen_weights_new,trial_data,bias_weights_new = self.chain1.calc_rosenbluth(UBase,retrace=False)
    if trial_data['abort']:
      mc_move_data['string'] = 'bad_trial_move_new'
      accept = False
      return accept,mc_move_data

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    rosen_weights_old,old_trial_data,bias_weights_old = self.chain1.calc_rosenbluth(UBase,retrace=True)
    if old_trial_data['abort']:
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
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))*np.product(bias_weights_old)
    Wold = np.product(np.sum(rosen_weights_old,axis=1))*np.product(bias_weights_new)
    mc_move_data['Wnew'] = Wnew
    mc_move_data['Wold'] = Wold
    mc_move_data['k'] = self.chain1.regrowth_length

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
      mc_move_data['U'] = trial_data['U']
      for local_index,sys_index in enumerate(self.chain1.regrowth_indices):
        x = trial_data['x'][local_index]
        y = trial_data['y'][local_index]
        z = trial_data['z'][local_index]
        self.system.x[sys_index] = x
        self.system.y[sys_index] = y
        self.system.z[sys_index] = z

        imx = trial_data['imx'][local_index]
        imy = trial_data['imy'][local_index]
        imz = trial_data['imz'][local_index]
        self.system.imx[sys_index] = imx
        self.system.imy[sys_index] = imy
        self.system.imz[sys_index] = imz

        self.system.neighbor_list.update_bead(sys_index,x=x,y=y,z=z)

    mc_move_data['string'] = 'Wnew/Wold: {:3.2e}/{:3.2e}={:5.4f}'.format(Wnew,Wold,Wnew/Wold)
    return accept,mc_move_data









