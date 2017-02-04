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

class SeverSC_CBMCMove(MonteCarloMove):
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
  def __init__(self,regrowth_types,bias_pkl,num_trials=30,regrowth_min=4,regrowth_max=20,viz=None):
    super(SeverSC_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='SeverSC_CBMCMove'
    self.regrowth_types = regrowth_types

    with open(bias_pkl,'rb') as f:
      self.bias = cPickle.load(f)

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max

    # We'll pre-instatiate the rosen-chain so 
    self.rosen_chain = RosenbluthChain(num_trials,regrowth_min,regrowth_max,self.bias,viz=viz)
  def set_engine(self,engine):
    self.engine                  = engine
    self.system                  = engine.system
    self.rosen_chain.engine      = engine
    self.rosen_chain.system      = engine.system
  @MonteCarloMove.counter
  def attempt(self):
    mc_move_data = {'string':'SEVR::'}

    abort,mol,counts = self.pick_mol()

    if abort:
      accept=False
      mc_move_data['string'] += 'not_enough_loops_or_tails'
      return accept,mc_move_data

    self.rosen_chain.chain_type = 'sever'
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
    abort,rosen_weights_new,bias_weights_new,Jnew = self.rosen_chain.calc_rosenbluth(UBase,retrace=False)
    if abort:
      mc_move_data['string'] += 'bad_trial_move_new'
      accept = False
      self.rosen_chain.reset()
      return accept,mc_move_data

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    abort,rosen_weights_old,bias_weights_old,Jold = self.rosen_chain.calc_rosenbluth(UBase,retrace=True)
    if abort:
      mc_move_data['string'] += 'bad_trial_move_old'
      accept = False
      self.rosen_chain.reset()
      return accept,mc_move_data

    #################################
    ## CALCULATE FULL ROSEN FACTOR ##
    #################################
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))*np.product(bias_weights_old)*Jnew
    Wold = np.product(np.sum(rosen_weights_old,axis=1))*np.product(bias_weights_new)*Jold
    mc_move_data['Wnew'] = Wnew
    mc_move_data['Wold'] = Wold
    mc_move_data['string'] += 'Wnew/Wold: {:3.2e}/{:3.2e}={:5.4f}'.format(Wnew,Wold,Wnew/Wold)

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

    # if accept:
    #   print accept,mc_move_data['string']
    #   self.engine.viz.clear()
    #   self.engine.viz.draw_system(bonds=True)
    #   self.engine.viz.show()
    #   ist()

    self.rosen_chain.reset()
    return accept,mc_move_data
  def pick_mol(self):
    tail_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tail']
    tie_list  = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tie']
    loop_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='loop']

    if len(tie_list)==0 and len(loop_list)==0:
      abort = True
      counts = {'tail':0,'tie':0,'loop':0}
      return abort,None,counts

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tie_list for i in mol.indices]
    index_list_loop = [i for mol in loop_list for i in mol.indices]
    counts = {'tail':len(index_list_tail),'tie':len(index_list_tie),'loop':len(index_list_loop)}

    index_list = index_list_loop + index_list_tie

    index = choice(index_list)
    mol = self.system.molecule_map[index]
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

    self.rosen_chain.set_indices(indices,internal_growth=True,random_inversion=False)
    self.rosen_chain.build_arrays()
    self.rosen_chain.copy_retrace_to_regrowth()

    chain1_end_sys   = self.rosen_chain.indices[-1]
    chain2_start_sys = self.rosen_chain.high_indices[0]
    l1 = self.rosen_chain.length

    self.rosen_chain.regrowth_beacons = [None for i in range(l1)] 
    self.rosen_chain.regrowth_anchors[-1] = (self.rosen_chain.regrowth_anchors[-1][1],)

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['removed'].append([chain1_end_sys,chain2_start_sys])

    mod_mol1 = {}
    mod_mol1['molecule']  = mol1
    mod_mol1['indices']   = self.rosen_chain.low_indices + self.rosen_chain.indices

    mol2 = ChainSegment()
    mol2.indices = self.rosen_chain.high_indices

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[mod_mol1],'added':[mol2],'removed':[]}
    # print "CHAIN1:"
    # print mol1.indices
    # print "GROW1:"
    # print self.rosen_chain.indices
    # print "FULL BONDS"
    # print self.rosen_chain.regrowth_bonds
    # print "BOND DELTA"
    # print bonds
    # print "ACC PKG"
    # print self.rosen_chain.acceptance_package
    # ist()
