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
  def __init__(self,beacon_file,num_trials=30,regrowth_min=4,regrowth_max=20,viz=None,bias=None):
    super(SeverSC_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='SeverSC_CBMCMove'
    self.bias = bias

    with open(beacon_file,'rb') as f:
      beacon = cPickle.load(f)

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max

    # We'll pre-instatiate the rosen-chain so 
    self.rosen_chain = RosenbluthChain(num_trials,regrowth_min,regrowth_max,beacon,viz=viz)
  def set_engine(self,engine):
    self.engine                  = engine
    self.system                  = engine.system
    self.rosen_chain.engine      = engine
    self.rosen_chain.system      = engine.system
  def _attempt(self):
    self.reset('SEVR::')

    abort,mol,counts = self.pick_mol()

    if abort:
      self.technical_abort = True
      self.accept          = False
      self.string += 'not_enough_loops_or_tails'
      return 

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


    ###################################
    ## CALCULATE TOPOLOGICAL CHANGES ##
    ###################################
    loopsNew = counts['loops'] + self.rosen_chain.acceptance_package['delta_topology']['loops']
    tiesNew  = counts['ties']  + self.rosen_chain.acceptance_package['delta_topology']['ties']
    tailsNew = counts['tails'] + self.rosen_chain.acceptance_package['delta_topology']['tails']
    loopsOld = counts['loops'] 
    tiesOld  = counts['ties']  
    tailsOld = counts['tails'] 

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
    Cnew = 1
    Bnew = self.bias.factor(**self.rosen_chain.acceptance_package['delta_topology'])
    Wnew = Wnew * Gnew * Jnew * Cnew * Bnew

    Wold = np.product(np.sum(rosen_weights_old,axis=1))
    Gold = np.product(beacon_weights_new) #new is correct
    Jold = Jold
    Cold = 1
    Bold = 1#self.bias.factor(ties=tiesOld,tails=tailsOld,loops=loopsOld)
    Wold = Wold * Gold * Jold * Cold * Bold

    # print'===================================================='
    # print 'SEVR::Wnew | Wnew_us',Wnew,'|',Wnew/Bnew
    # print 'SEVR::Wold | Wold_us',Wold,'|',Wold/Bold
    # print 'SEVR::WRatio',Wnew/Wold, '|',Wnew*Bold/Wold/Bnew
    # print 'SEVR::Bnew,Bold',Bnew,Bold,Bnew/Bold
    # print 'SEVR::OLD::ties,loops,tails',tiesOld,tailsOld,loopsOld
    # print 'SEVR::NEW::ties,loops,tails',tiesNew,tailsNew,loopsNew
    # ist()

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

    if len(tie_list)==0 and len(loop_list)==0:
      abort = True
      return abort,None,counts

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tie_list for i in mol.indices]
    index_list_loop = [i for mol in loop_list for i in mol.indices]
    for name,l in zip(['tail_beads','tie_beads','loop_beads'],[index_list_tail,index_list_tie,index_list_loop]):
      counts[name] = len(l)

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

    if mol1.properties['topology'] == 'tie':
      tie_change  = -1
      tail_change = +2
      loop_change =  0
    elif mol1.properties['topology'] == 'loop':
      tie_change  =  0
      tail_change = +2
      loop_change = -1
    else:
      raise ValueError('Somehow, we have a tail chosen for a SeverSC move...')
    self.rosen_chain.acceptance_package['delta_topology'] = {'ties':tie_change, 'loops':loop_change,'tails':tail_change}
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
