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
  def __init__(self,regrowth_types,bias_pkl,num_trials=30,regrowth_min=4,regrowth_max=10,viz=None):
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
    self.rosen_chain.engine      = engine
    self.rosen_chain.system      = engine.system
  def _attempt(self):
    self.reset('SPLC::')

    abort,mol1,mol2,counts = self.pick_mol()

    if abort:
      self.technical_abort = True
      self.accept=False
      if mol1 is None:
        self.string += 'not_enough_tails'
      elif mol2 is None:
        self.string += 'chain1_too_short'
      else:
        self.string += 'chain2_too_short'
      return

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
    abort,rosen_weights_new,bias_weights_new,Jnew = self.rosen_chain.calc_rosenbluth(UBase,retrace=False)
    if abort:
      self.string += 'bad_trial_move_new'
      self.accept = False
      self.rosen_chain.reset()
      return 

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    abort,rosen_weights_old,bias_weights_old,Jold = self.rosen_chain.calc_rosenbluth(UBase,retrace=True)
    if abort:
      self.string += 'bad_trial_move_old'
      self.accept = False
      self.rosen_chain.reset()
      return

    #################################
    ## CALCULATE FULL ROSEN FACTOR ##
    #################################
    count_new = float(counts['num_close'])
    count_old = float(counts['ties'] + counts['loops'] + 1)
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))*np.product(bias_weights_old)*Jnew# *1/count_new
    Wold = np.product(np.sum(rosen_weights_old,axis=1))*np.product(bias_weights_new)*Jold# *1/count_old

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
    #   print accept,mc_move_data['string']
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

    if len(tail_list)<2:
      abort = True
      return abort,None,None,counts

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tie_list for i in mol.indices]
    index_list_loop = [i for mol in loop_list for i in mol.indices]
    for name,l in zip(['tail_beads','tie_beads','loop_beads'],[index_list_tail,index_list_tie,index_list_loop]):
      counts[name] = len(l)

    index_list_tail_ends = []
    for mol in tail_list: 
      for i in mol.properties['chain_ends']:
        if i not in mol.properties['connected_to']:
          index_list_tail_ends.append(i)

    index_list = index_list_tail_ends[:]

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
    close_mask = (index_list_close!=-1)
    index_list_close = index_list_close[close_mask]
    index2 = index_list[choice(index_list_close)]
    mol2 = self.system.molecule_map[index2]
    counts['num_close'] = len(index_list_close)

    # print mol1,mol2,index1,index2
    # print 'MOL1/MOL2 INDICES'
    # print mol1.indices
    # print mol2.indices
    # print 'MOL1/MOL2 PROPERTIES'
    # print mol1.properties
    # print mol2.properties


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

    self.rosen_chain.set_indices(indices1,internal_growth=False,random_inversion=False)#,full_chain=True)
    self.rosen_chain.build_arrays()
    self.rosen_chain.copy_retrace_to_regrowth()

    chain1_end_sys   = self.rosen_chain.indices[-1]
    chain1_end_new   = self.rosen_chain.trial_indices[-1]
    chain2_start_sys = indices2[0]

    l1 = self.rosen_chain.length
    self.rosen_chain.regrowth_beacons = [(chain2_start_sys,l1-i) for i in range(l1)] 
    self.rosen_chain.regrowth_anchors[-1] += (chain2_start_sys,)


    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([chain1_end_sys,chain2_start_sys])

    mod_mol1 = {}
    mod_mol1['molecule']     = mol1
    mod_mol1['indices'] = indices1 + indices2

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[mod_mol1],'added':[],'removed':[mol2]}
