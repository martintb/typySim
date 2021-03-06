import cPickle
import numpy as np
from numpy.random import randint,random,choice 
import ipdb as pdb; ist = pdb.set_trace

from typySim import molecule
from typySim.geometry import linalg
from typySim.engine.RosenbluthChain import RosenbluthChain
from typySim.engine.MonteCarloMove import MonteCarloMove
from typySim.core.cy.cyutil import n_closest

class SwapSC_CBMCMove(MonteCarloMove):
  '''
  Definitions
  -----------
  local index :  [0,len(regrowth_mol.indices))
    Integer ranging between 0 and the length of the chain being regrown. Used
    to access the indices and data arrays for the regrowing chain.

  sys index : [0,self.system.nbeads)
    Global index of a bead in the system arrays (e.g. self.system.x ). This value ranges from .  
  new index : [self.system.nbeads,self.system.nbeads+len(regrowth_mol.indices))
    Global index of the trial_move bead. This value ranges from self.system to
    the number of beads being regrown

  XXX_index_trial : [0,self.num_trials)
    Represents the index of one of the possibi

  '''
  def __init__(self,beacon_file,num_trials=30,regrowth_min=4,regrowth_max=10,viz=None,bias=None):
    super(SwapSC_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='SwapSC_CBMCMove'
    self.bias = bias

    with open(beacon_file,'rb') as f:
      beacon = cPickle.load(f)

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max

    # We'll pre-instatiate the rosen-chain so 
    self.rosen_subchain1 = RosenbluthChain(num_trials,regrowth_min,regrowth_max,beacon,viz=None)
    self.rosen_subchain2 = RosenbluthChain(num_trials,regrowth_min,regrowth_max,beacon,viz=None)
    self.rosen_chain = RosenbluthChain(num_trials,regrowth_min,regrowth_max,beacon,viz=viz)
  def set_engine(self,engine):
    self.engine                  = engine
    self.system                  = engine.system
    self.rosen_chain.engine      = engine
    self.rosen_chain.system      = engine.system
    self.rosen_subchain1.engine  = engine
    self.rosen_subchain1.system  = engine.system
    self.rosen_subchain2.engine  = engine
    self.rosen_subchain2.system  = engine.system
  def _attempt(self):
    self.reset('SWAP::')

    # abort,mol1,mol2,counts = self.pick_mol()
    # if abort:
    #   self.technical_abort = True
    #   self.accept=False
    #   if mol1 is None: 
    #     self.string += 'not_enough_tails'
    #   elif mol2 is None:
    #     self.string += 'chain1_too_short'
    #   else:
    #     self.string += 'chain2_too_short'
    #   return
    # self.rosen_chain.chain_type = 'splice'
    # self.build_rosen_chain(mol1,mol2)

    result = self.pick_mol2()
    if result[0] == True:
      self.technical_abort = True
      self.accept=False
      self.string += result[1]
      return
    _,mol1,mol2,counts,index1,index2 = result
    self.build_rosen_chain(mol1,mol2)
    # self.build_rosen_chain2(mol1,mol2,index1,index2)

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
      self.rosen_subchain1.reset()
      self.rosen_subchain2.reset()
      return

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    abort,rosen_weights_old,beacon_weights_old,Jold = self.rosen_chain.calc_rosenbluth(UBase,retrace=True)
    if abort:
      self.string += 'bad_trial_move_old'
      self.accept = False
      self.rosen_chain.reset()
      self.rosen_subchain1.reset()
      self.rosen_subchain2.reset()
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
    Cnew = 1
    Bnew = self.bias.factor(**self.rosen_chain.acceptance_package['delta_topology'])
    Wnew = Wnew * Gnew * Jnew * Cnew * Bnew

    Wold = np.product(np.sum(rosen_weights_old,axis=1))
    Gold = np.product(beacon_weights_new) #new is correct
    Jold = Jold
    Cold = 1
    Bold = 1 #self.bias.factor(ties=tiesOld,tails=tailsOld,loops=loopsOld)
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
      # print self.string
      # self.engine.viz.clear()
      # self.engine.viz.draw_system(bonds=True)
      # self.engine.viz.show()
      # ist()

      self.Unew = self.rosen_chain.acceptance_package['U']
      self.rosen_chain.apply_acceptance_package()

    self.string += 'Wnew/Wold: {:3.2e}/{:3.2e}={:5.4f}'.format(Wnew,Wold,Wnew/Wold)
    self.rosen_chain.reset()
    self.rosen_subchain1.reset()
    self.rosen_subchain2.reset()
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

    

    # index_list = index_list_tail + index_list_tie + index_list_loop
    mol_list = tail_list + tie_list + loop_list

    mol1 = choice(mol_list)
    mol_list.remove(mol1)
    
    if mol1.size < self.regrowth_min+2:
      abort = True
      return abort,mol1,None,counts

    x1,y1,z1 = mol1.properties['center_of_mass']
    r2 = [mol.properties['center_of_mass'] for mol in mol_list]
    x2Array,y2Array,z2Array = np.array(r2).T

    index_list_close,dist_list_close = n_closest(10,x1,y1,z1,x2Array,y2Array,z2Array,self.system.box)
    close_mask = (index_list_close!=-1)
    index_list_close = index_list_close[close_mask]
    mol2 = mol_list[choice(index_list_close)]
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
  def pick_mol2(self):
    counts = {}

    tail_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tail']
    tie_list  = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tie']
    loop_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='loop']
    total_count = 0
    for name,l in zip(['tails','ties','loops'],[tail_list,tie_list,loop_list]):
      counts[name] = len(l)
      total_count += len(l)

    if total_count<2:
      abort = True
      return abort,'not_enough_chains'

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tie_list for i in mol.indices]
    index_list_loop = [i for mol in loop_list for i in mol.indices]
    for name,l in zip(['tail_beads','tie_beads','loop_beads'],[index_list_tail,index_list_tie,index_list_loop]):
      counts[name] = len(l)

    index_list = index_list_tail + index_list_tie + index_list_loop

    index1 = choice(index_list)
    mol1 = self.system.molecule_map[index1]
    if mol1.size < self.regrowth_min+2:
      abort = True
      return abort,'chain1_too_short'

    index2_list = self.system.box.neighbor_list.get_neighbors_by_index(index1)
    index2_list = [i for i in index2_list if (i not in mol1.indices) and (self.system.molecule_map[i].name=='ChainSegment')]
    if len(index2_list)==0:
      abort = True
      return abort,'no_nearby_beads'

    x1 = self.system.x[index1]
    y1 = self.system.y[index1]
    z1 = self.system.z[index1]
    x2Array = np.array([self.system.x[i] for i in index2_list])
    y2Array = np.array([self.system.y[i] for i in index2_list])
    z2Array = np.array([self.system.z[i] for i in index2_list])

    index_list_close,dist_list_close = n_closest(10,x1,y1,z1,x2Array,y2Array,z2Array,self.system.box)
    close_mask = (index_list_close!=-1)
    index_list_close = index_list_close[close_mask]
    counts['num_close'] = len(index_list_close)

    index2 = index2_list[choice(index_list_close)]
    mol2 = self.system.molecule_map[index2]

    # print mol1,mol2,index1,index2
    # print 'MOL1/MOL2 INDICES'
    # print mol1.indices
    # print mol2.indices
    # print 'MOL1/MOL2 PROPERTIES'
    # print mol1.properties
    # print mol2.properties


    if mol2.size < self.regrowth_min+2:
      abort = True
      return abort,'chain2_too_short'

    abort = False
    return abort,mol1,mol2,counts,index1,index2
  def build_rosen_chain(self,mol1,mol2):
    '''

    '''


    indices1 = mol1.indices
    if mol1.properties['topology'] == 'tail':
      # We want to be growing in the direction away from the crystallite
      if indices1[-1] in mol1.properties['connected_to']:
        indices1 = indices1[::-1]

    indices2 = mol2.indices
    if mol2.properties['topology'] == 'tail':
      # We want to be growing in the direction away from the crystallite
      if indices2[-1] in mol2.properties['connected_to']:
        indices2 = indices2[::-1]

    self.rosen_subchain1.set_indices(indices1,internal_growth=True,random_inversion=False)#,full_chain=True)
    self.rosen_subchain2.set_indices(indices2,internal_growth=True,random_inversion=False)#,full_chain=True)
    self.rosen_subchain1.build_arrays()
    self.rosen_subchain2.build_arrays(index_shift = self.rosen_subchain1.length)
    self.rosen_subchain1.copy_retrace_to_regrowth()
    self.rosen_subchain2.copy_retrace_to_regrowth()


    new_beacon1          = self.rosen_subchain2.retrace_beacons[0][0]
    new_beacon2          = self.rosen_subchain1.retrace_beacons[0][0]
    chain1_end_sys       = self.rosen_subchain1.indices[-1]
    chain2_end_sys       = self.rosen_subchain2.indices[-1]
    chain1_last_anchor   = self.rosen_subchain1.trial_indices[-2]
    chain2_last_anchor   = self.rosen_subchain2.trial_indices[-2]
    chain1_end_trial     = self.rosen_subchain1.trial_indices[-1]
    chain2_end_trial     = self.rosen_subchain2.trial_indices[-1]

    l1 = self.rosen_subchain1.length
    l2 = self.rosen_subchain2.length
    self.rosen_subchain1.regrowth_beacons = [(new_beacon1,l1-i) for i in range(l1)] 
    self.rosen_subchain2.regrowth_beacons = [(new_beacon2,l2-i) for i in range(l2)] 
    self.rosen_subchain1.regrowth_anchors[-1] = (new_beacon1,chain1_last_anchor)
    self.rosen_subchain2.regrowth_anchors[-1] = (new_beacon2,chain2_last_anchor)

    self.rosen_chain.combine(self.rosen_subchain1,self.rosen_subchain2)

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([chain1_end_sys,new_beacon1])
    bonds['added'].append([chain2_end_sys,new_beacon2])
    bonds['removed'].append([chain1_end_sys,new_beacon2])
    bonds['removed'].append([chain2_end_sys,new_beacon1])

    mod_mol1 = {}
    mod_mol1['molecule']     = mol1
    mod_mol1['indices'] = self.rosen_subchain1.low_indices + self.rosen_subchain1.indices + self.rosen_subchain2.high_indices

    mod_mol2 = {}
    mod_mol2['molecule']     = mol2
    mod_mol2['indices'] = self.rosen_subchain2.low_indices + self.rosen_subchain2.indices + self.rosen_subchain1.high_indices

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[mod_mol1,mod_mol2],'added':[],'removed':[]}

    tie_change  = 0
    tail_change = 0
    loop_change = 0
    for mol in [mol1,mol2]:
      if mol.properties['topology'] == 'tie':
        tie_change   -= 1
      elif mol.properties['topology'] == 'loop':
        loop_change  -= 1
      elif mol.properties['topology'] == 'tail':
        tail_change  -= 1


    mol1_connect =     [mol1.properties['connected_to'].get(mod_mol1['indices'][0] ,None)]
    mol1_connect.append(mol2.properties['connected_to'].get(mod_mol1['indices'][-1],None))
    mol2_connect =     [mol2.properties['connected_to'].get(mod_mol2['indices'][0] ,None)]
    mol2_connect.append(mol1.properties['connected_to'].get(mod_mol2['indices'][-1],None))

    for cmol1,cmol2 in [mol1_connect,mol2_connect]:
      if (cmol1 is None) or (cmol2 is None):
        tail_change += 1
      elif cmol1['molecule'] is cmol2['molecule']:
        loop_change += 1
      elif cmol1['molecule'] is not cmol2['molecule']:
        tie_change  += 1
      elif (cmol1 is None) and (cmol2 is None):
        raise ValueError('Something is wrong with SwapSC topology shift math')
    
    self.rosen_chain.acceptance_package['delta_topology'] = {'ties':tie_change, 'loops':loop_change,'tails':tail_change}
  def build_rosen_chain2(self,mol1,mol2,index1,index2):
    '''

    '''


    indices1 = mol1.indices
    if mol1.properties['topology'] == 'tail':
      # We want to be growing in the direction away from the crystallite
      if indices1[-1] in mol1.properties['connected_to']:
        indices1 = indices1[::-1]

    indices2 = mol2.indices
    if mol2.properties['topology'] == 'tail':
      # We want to be growing in the direction away from the crystallite
      if indices2[-1] in mol2.properties['connected_to']:
        indices2 = indices2[::-1]

    self.rosen_subchain1.set_indices_around(indices1,index1,random_inversion=False)#,full_chain=True)
    self.rosen_subchain2.set_indices_around(indices2,index2,random_inversion=False)#,full_chain=True)
    self.rosen_subchain1.build_arrays()
    self.rosen_subchain2.build_arrays(index_shift = self.rosen_subchain1.length)
    self.rosen_subchain1.copy_retrace_to_regrowth()
    self.rosen_subchain2.copy_retrace_to_regrowth()


    new_beacon1          = self.rosen_subchain2.retrace_beacons[0][0]
    new_beacon2          = self.rosen_subchain1.retrace_beacons[0][0]
    chain1_end_sys       = self.rosen_subchain1.indices[-1]
    chain2_end_sys       = self.rosen_subchain2.indices[-1]
    chain1_last_anchor   = self.rosen_subchain1.trial_indices[-2]
    chain2_last_anchor   = self.rosen_subchain2.trial_indices[-2]
    chain1_end_trial     = self.rosen_subchain1.trial_indices[-1]
    chain2_end_trial     = self.rosen_subchain2.trial_indices[-1]

    l1 = self.rosen_subchain1.length
    l2 = self.rosen_subchain2.length
    self.rosen_subchain1.regrowth_beacons = [(new_beacon1,l1-i) for i in range(l1)] 
    self.rosen_subchain2.regrowth_beacons = [(new_beacon2,l2-i) for i in range(l2)] 
    self.rosen_subchain1.regrowth_anchors[-1] = (new_beacon1,chain1_last_anchor)
    self.rosen_subchain2.regrowth_anchors[-1] = (new_beacon2,chain2_last_anchor)

    self.rosen_chain.combine(self.rosen_subchain1,self.rosen_subchain2)

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([chain1_end_sys,new_beacon1])
    bonds['added'].append([chain2_end_sys,new_beacon2])
    bonds['removed'].append([chain1_end_sys,new_beacon2])
    bonds['removed'].append([chain2_end_sys,new_beacon1])

    mod_mol1 = {}
    mod_mol1['molecule']     = mol1
    mod_mol1['indices'] = self.rosen_subchain1.low_indices + self.rosen_subchain1.indices + self.rosen_subchain2.high_indices

    mod_mol2 = {}
    mod_mol2['molecule']     = mol2
    mod_mol2['indices'] = self.rosen_subchain2.low_indices + self.rosen_subchain2.indices + self.rosen_subchain1.high_indices

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'modded':[mod_mol1,mod_mol2],'added':[],'removed':[]}

    tie_change  = 0
    tail_change = 0
    loop_change = 0
    for mol in [mol1,mol2]:
      if mol.properties['topology'] == 'tie':
        tie_change   -= 1
      elif mol.properties['topology'] == 'loop':
        loop_change  -= 1
      elif mol.properties['topology'] == 'tail':
        tail_change  -= 1


    mol1_connect =     [mol1.properties['connected_to'].get(mod_mol1['indices'][0] ,None)]
    mol1_connect.append(mol2.properties['connected_to'].get(mod_mol1['indices'][-1],None))
    mol2_connect =     [mol2.properties['connected_to'].get(mod_mol2['indices'][0] ,None)]
    mol2_connect.append(mol1.properties['connected_to'].get(mod_mol2['indices'][-1],None))

    for cmol1,cmol2 in [mol1_connect,mol2_connect]:
      if (cmol1 is None) or (cmol2 is None):
        tail_change += 1
      elif cmol1['molecule'] is cmol2['molecule']:
        loop_change += 1
      elif cmol1['molecule'] is not cmol2['molecule']:
        tie_change  += 1
      elif (cmol1 is None) and (cmol2 is None):
        raise ValueError('Something is wrong with SwapSC topology shift math')
    
    self.rosen_chain.acceptance_package['delta_topology'] = {'ties':tie_change, 'loops':loop_change,'tails':tail_change}
