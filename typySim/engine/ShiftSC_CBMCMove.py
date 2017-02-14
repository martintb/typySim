import cPickle
import numpy as np
from numpy.random import randint,random,choice
import ipdb as pdb; ist = pdb.set_trace

from typySim import molecule
from typySim.geometry import linalg
from typySim.engine.RosenbluthChain import RosenbluthChain
from typySim.engine.MonteCarloMove import MonteCarloMove
from typySim.molecule.ChainSegment import ChainSegment

from typySim.core.cy.cyutil import n_closest

class ShiftSC_CBMCMove(MonteCarloMove):
  '''
  The ShiftSC move is a standard CBMC Regrowth move that either regrows the entire chain
  or an internal segment of the chain.


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
  def __init__(self,beacon_file,num_trials=30,regrowth_min=4,regrowth_max=20,viz=None,bias=None):
    super(ShiftSC_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='ShiftSC_CBMCMove'
    self.bias = None

    with open(beacon_file,'rb') as f:
      self.beacon = cPickle.load(f)

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max

    # We'll pre-instatiate the rosen-chain so 
    self.rosen_chain = RosenbluthChain(num_trials,regrowth_min,regrowth_max,self.beacon,viz=viz)
  def set_engine(self,engine):
    self.engine                  = engine
    self.system                  = engine.system
    self.rosen_chain.engine      = engine
    self.rosen_chain.system      = engine.system
  def _attempt(self):
    self.reset('SHFT::')

    abort,mol,counts = self.pick_mol()

    if abort:
      self.technical_abort = True
      self.accept          = False
      if mol is None:
        self.string += 'no_sc_chains'
      else:
        self.string += 'chain_too_short'
      return

    self.rosen_chain.chain_type = 'shift'
    self.build_rosen_chain(mol)

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
    abort,rosen_weights_new,beacon_weights_new,Jnew = self.rosen_chain.calc_rosenbluth(UBase,retrace=False)
    if abort:
      self.string += 'bad_trial_move_new'
      self.accept = False
      self.rosen_chain.reset()
      return 

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    abort,rosen_weights_old,beacon_weights_old,Jold = self.rosen_chain.calc_rosenbluth(UBase,retrace=True)
    if abort:
      self.string += 'bad_trial_move_old'
      self.accept = False
      self.rosen_chain.reset()
      return 

    #################################
    ## CALCULATE FULL ROSEN FACTOR ##
    #################################
    # W: base rosenbluth factor of the configuration
    # G: beacon factor used in fixed endpoint_CBMC
    # J: constraint factor used for last step of fixed endpoint CBMC
    # C: count bias based on how we selected the chain to peturb
    # B: artificial bias used to control the 
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))
    Gnew = np.product(beacon_weights_old) #old is correct
    Jnew = Jnew
    Cnew = 1 # chance of choosing chain is identical in both directions
    Bnew = 1 # no topology change so no topology change bias
    Wnew = Wnew * Gnew * Jnew * Cnew * Bnew

    Wold = np.product(np.sum(rosen_weights_old,axis=1))
    Gold = np.product(beacon_weights_new) #new is correct
    Jold = Jold
    Cold = 1 # chance of choosing chain is identical in both directions
    Bold = 1 # no topology change so no topology change bias
    Wold = Wold * Gold * Jold * Cold * Bold

    #######################
    ## ACCEPT OR REJECT? ##
    #######################
    if Wnew>Wold:
      self.accept = True
    else:
      ranf = np.random.rand()
      if ranf<(Wnew/Wold):
        self.accept=True
      else:
        self.accept=False

    if self.accept:
      self.Unew = self.rosen_chain.acceptance_package['U']
      self.rosen_chain.apply_acceptance_package()

    # if accept:
    #   print self.string
    #   self.engine.viz.clear()
    #   self.engine.viz.draw_system(bonds=True)
    #   self.engine.viz.show()
    #   ist()

    self.string += 'Wnew/Wold: {:3.2e}/{:3.2e}={:5.4f}'.format(Wnew,Wold,Wnew/Wold)
    self.rosen_chain.reset()
    return 
  def pick_mol(self):
    counts = {}

    tail_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tail']
    tie_list  = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tie']
    loop_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='loop']
    for name,l in zip(['tails','ties','loops'],[tail_list,tie_list,loop_list]):
      counts[name] = len(l)

    if len(tie_list)==0 and len(loop_list)==0 and len(tail_list)==0:
      abort = True
      return abort,None,counts

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tie_list for i in mol.indices]
    index_list_loop = [i for mol in loop_list for i in mol.indices]
    for name,l in zip(['tail_beads','tie_beads','loop_beads'],[index_list_tail,index_list_tie,index_list_loop]):
      counts[name] = len(l)

    index_list = index_list_loop + index_list_tie + index_list_tail

    index = choice(index_list)
    mol = self.system.molecule_map[index]

    if mol.size < self.regrowth_min+2:
      abort = True
      return abort,mol,counts

    abort = False
    return abort,mol,counts
  def build_rosen_chain(self,mol1):
    '''
      >> Pick a loop or tails to break up
      S1--00--NR--RR--RR--RR--RR--RR--RR--RR
                                           |
      S2--10--11--12--13--14--15--16--17--RR

      >> Pick portion of chain to regrow (RR = beads being regrown)
      >> Regrow using 00 as the first step anchor, and 17 as the beacon for all steps.
      S1--00--RR--RR--RR--RR--RR--RR--RR--RR--RR
      S2--10--11--12--13--14--15--16--17

    '''
    indices = mol1.indices

    if mol1.properties['topology'] == 'tail':
      # We want to be growing in the direction away from the crystallite
      if indices[-1] in mol1.properties['connected_to']:
        indices = indices[::-1]

      self.rosen_chain.set_indices(indices,internal_growth=False,random_inversion=False)
    else:
      self.rosen_chain.set_indices(indices,internal_growth=True)

    self.rosen_chain.build_arrays()
    self.rosen_chain.copy_retrace_to_regrowth()
    self.rosen_chain.acceptance_package['delta_topology'] = {'ties':0, 'loops':0,'tails':0}
    self.rosen_chain.acceptance_package['bonds'] = {'added':[], 'removed':[]}
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[],'added':[],'removed':[]}
