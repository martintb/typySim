import cPickle
import numpy as np
from numpy.random import randint,random,choice
import ipdb as pdb; ist = pdb.set_trace

from typySim import molecule
from typySim.geometry import linalg
from typySim.engine.RosenbluthChain import RosenbluthChain
from typySim.engine.MonteCarloMove import MonteCarloMove
from typySim.core.cy.cyutil import n_closest

class SpliceSC_CBMCMove(MonteCarloMove):
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
    super(SpliceSC_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='SpliceSC_CBMCMove'
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
    self.rosen_subchain1.engine  = engine
    self.rosen_subchain1.system  = engine.system
    self.rosen_subchain2.engine  = engine
    self.rosen_subchain2.system  = engine.system
    self.rosen_chain.engine      = engine
    self.rosen_chain.system      = engine.system
  @MonteCarloMove.counter
  def attempt(self):
    mc_move_data = {'string':'SPLC::'}

    abort,mol1,mol2,counts = self.pick_mol()

    if abort:
      accept=False
      if mol1 is None:
        mc_move_data['string'] += 'not_enough_tails'
        return accept,mc_move_data
      elif mol2 is None:
        mc_move_data['string'] += 'chain1_too_short'
        return accept,mc_move_data
      else:
        mc_move_data['string'] = 'chain2_too_short'
        return accept,mc_move_data

    self.rosen_chain.chain_type = 'splice'
    self.build_rosen_chain(mol1,mol2)

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
    abort,rosen_weights_new,bias_weights_new = self.rosen_chain.calc_rosenbluth(UBase,retrace=False)
    if abort:
      mc_move_data['string'] += 'bad_trial_move_new'
      accept = False
      self.rosen_subchain1.reset()
      self.rosen_subchain2.reset()
      self.rosen_chain.reset()
      return accept,mc_move_data

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    abort,rosen_weights_old,bias_weights_old = self.rosen_chain.calc_rosenbluth(UBase,retrace=True)
    if abort:
      mc_move_data['string'] += 'bad_trial_move_old'
      accept = False
      self.rosen_subchain1.reset()
      self.rosen_subchain2.reset()
      self.rosen_chain.reset()
      return accept,mc_move_data

    #################################
    ## CALCULATE FULL ROSEN FACTOR ##
    #################################
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))*np.product(bias_weights_old)
    Wold = np.product(np.sum(rosen_weights_old,axis=1))*np.product(bias_weights_new)
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

    if accept:
      print accept,mc_move_data['string']
      self.engine.viz.clear()
      self.engine.viz.draw_system(bonds=True)
      self.engine.viz.show()
      ist()

    self.rosen_subchain1.reset()
    self.rosen_subchain2.reset()
    self.rosen_chain.reset()
    return accept,mc_move_data
  def pick_mol(self):
    tail_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tail']
    tie_list  = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tie']
    loop_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='loop']

    if len(tail_list)<2:
      abort = True
      counts = {'tail':0,'tie':0,'loop':0}
      return abort,None,None,counts

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tail_list for i in mol.indices]
    index_list_loop = [i for mol in tail_list for i in mol.indices]
    counts = {'tail':len(index_list_tail),'tie':len(index_list_tie),'loop':len(index_list_loop)}

    index_list = index_list_tail[:]

    index1 = choice(index_list)
    mol1 = self.system.molecule_map[index1]
    index_list = [i for i in index_list if i not in mol1.indices]

    if mol1.size < self.regrowth_min+2:
      abort = True
      return abort,mol1,None,counts

    x1 = self.system.x[index1]
    y1 = self.system.y[index1]
    z1 = self.system.z[index1]
    x2Array = self.system.x[index_list]
    y2Array = self.system.y[index_list]
    z2Array = self.system.z[index_list]

    index_list_close,dist_list_close = n_closest(10,x1,y1,z1,x2Array,y2Array,z2Array,self.system.box)
    index2 = index_list[choice(index_list_close)]
    mol2 = self.system.molecule_map[index2]

    if mol2.size < self.regrowth_min+2:
      abort = True
      return abort,mol1,mol2,counts

    abort = False
    return abort,mol1,mol2,counts
  def build_rosen_chain(self,mol1,mol2):
    '''
    There are several approaches you could use to set up the splice move. This is 
    one approach that seems to be reasonably successful. 

      >> Pick two tails to regrow
      S1--00--01--02--03--04--05--06--07--08--09
      S2--10--11--12--13--14--15--16--17

      >> Pick longer chain to regrow (RR = beads being regrown)
      S1--00--RR--RR--RR--RR--RR--RR--RR--RR--RR
      S2--10--11--12--13--14--15--16--17

      >> Regrow using 00 as the first step anchor, and 17 as the beacon for all steps.
      S1--00--NR--RR--RR--RR--RR--RR--RR--RR
                                           |
      S2--10--11--12--13--14--15--16--17--RR

    '''
    if mol1.size<mol2.size:
      mol1,mol2 = mol2,mol1

    indices1 = mol1.indices
    indices2 = mol2.indices

    # We want to be growing in the direction away from the crystallite
    if indices1[-1] in mol1.properties['connected_to']:
      indices1 = indices1[::-1]

    if indices2[-1] not in mol2.properties['connected_to']:
      indices2 = indices2[::-1]

    self.rosen_chain.set_indices(indices1,internal_growth=False,random_inversion=False,full_chain=True)
    self.rosen_chain.build_arrays()
    self.rosen_chain.copy_retrace_to_regrowth()

    l1 = self.rosen_chain.length
    self.rosen_chain.regrowth_beacons = [(indices2[0],l1-i) for i in range(l1)] 

    chain1_end_sys   = self.rosen_chain.indices[-1]
    chain1_end_new   = self.rosen_chain.new_indices[-1]
    chain2_start_sys = indices2[0]

    self.rosen_chain.regrowth_bonds[-1].append([chain1_end_new,chain2_start_sys])

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([chain1_end_sys,chain2_start_sys])

    mod_mol1 = {}
    mod_mol1['molecule']     = mol1
    mod_mol1['indices'] = indices1 + indices2

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[mod_mol1],'added':[],'removed':[mol2]}
