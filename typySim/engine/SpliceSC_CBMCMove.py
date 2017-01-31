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

import copy



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
  def __init__(self,regrowth_types,bias_pkl,num_trials=4,regrowth_min = 10,regrowth_max=12,viz=None):
    super(SpliceSC_CBMCMove,self).__init__() #must call parent class' constructor
    self.name='SpliceSC_CBMCMove'
    self.regrowth_types = regrowth_types

    with open(bias_pkl,'rb') as f:
      self.bias = cPickle.load(f)

    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max

    # We'll pre-instatiate three rosen-chains that we can combine and split at 
    # will. 
    self.rosen_subchain1 = RosenbluthChain(num_trials,regrowth_min,regrowth_max,self.bias)
    self.rosen_subchain2 = RosenbluthChain(num_trials,regrowth_min,regrowth_max,self.bias)
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
    mc_move_data = {}


    abort,mol1,mol2,counts = self.pick_mol()
    topo1 = mol1.properties['topology']
    topo2 = mol2.properties['topology']

    if random()>0.5:
      self.rosen_chain.chain_type = 'regrow'
    elif random()>0.5:
      self.rosen_chain.chain_type = 'rebridge'
    elif topo1 in ['loop','tie']:
      self.rosen_chain.chain_type = 'sever'
    elif topo1=='tail' and topo2=='tail':
      self.rosen_chain.chain_type = 'splice'

    ####################################
    ## DETERMINE CHAINS BEING REGROWN ##
    ####################################
    # Choose a chain molecule and a starting point
    N1,chain1_dict = self.pick_chain()
    if chain1_dict['length']<self.regrowth_min+2:
      accept = False
      mc_move_data['string'] = 'chain1_is_too_short'
      return accept,mc_move_data

    sever = False
    if (not regrow) and (chain1_dict['topology'] in ['loop','tie']) and (random()>0.5):
      sever   = True

    if (not regrow) and (not sever):
      N2,chain2_dict = self.pick_chain(restrict_list=[chain1_dict['mol']])
      if chain2_dict['length']<self.regrowth_min+2:
        accept = False
        mc_move_data['string'] = 'chain2_is_too_short'
        return accept,mc_move_data
    else:
      N2 = 0
      chain2_dict = None

    if regrow:
      mc_move_data['string'] = 'RGRW'
      self.regrow(chain1_dict)
      self.rosen_chain.chain_type = 'regrow'
    elif sever:
      mc_move_data['string'] = 'SEVR'
      self.sever(chain1_dict)
      self.rosen_chain.chain_type = 'sever'
    elif (chain1_dict['topology'] == chain2_dict['topology'] == 'tail'):
      mc_move_data['string'] = 'SPCE'
      self.splice(chain1_dict,chain2_dict)
      self.rosen_chain.chain_type = 'splice'
    else:
      mc_move_data['string'] = 'RBGE'
      self.rebridge(chain1_dict,chain2_dict)
      self.rosen_chain.chain_type = 'rebridge'

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
      mc_move_data['string'] += '::bad_trial_move_new'
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
      mc_move_data['string'] += '::bad_trial_move_old'
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
      pkg = self.rosen_chain.acceptance_package
      mc_move_data['U'] = pkg['U']

      #------------------------#
      # UPDATE COORD AND IMAGE #
      #------------------------#
      for local_index,sys_index in enumerate(self.rosen_chain.indices):
        x = pkg['x'][local_index]
        y = pkg['y'][local_index]
        z = pkg['z'][local_index]
        self.system.x[sys_index] = x
        self.system.y[sys_index] = y
        self.system.z[sys_index] = z

        imx = pkg['imx'][local_index]
        imy = pkg['imy'][local_index]
        imz = pkg['imz'][local_index]
        self.system.imx[sys_index] = imx
        self.system.imy[sys_index] = imy
        self.system.imz[sys_index] = imz

        self.system.neighbor_list.update_bead(sys_index,x=x,y=y,z=z)

      #--------------#
      # UPDATE BONDS #
      #--------------#
      if 'bonds' in pkg:
        for beadi,beadj in pkg['bonds']['added']:
          self.system.bonds.add(beadi,beadj,0)
        for beadi,beadj in pkg['bonds']['removed']:
          self.system.bonds.remove(beadi,beadj,0)

      #------------------#
      # UPDATE MOLECULES #
      #------------------#
      if 'molecules' in pkg:
        for mol in pkg['molecules']['removed']:
          self.system.remove_molecule(molecule=mol,remove_beads=False)

        if 'mol1' in pkg['molecules']:
          for mol_str in ['mol1','mol2']:
            if mol_str not in pkg['molecules']:
              continue
            
            mol = pkg['molecules'][mol_str]['molecule']
            mol.indices = pkg['molecules'][mol_str]['indices']
            if mol in pkg['molecules']['added']:
              self.system.add_molecule(mol)
            else:
              for index in mol.indices:
                self.system.molecule_map[index] = mol

            mol.reset()
            mol.update_properties() #update topology
            if mol.properties['topology'] == 'tail':
              mol.types[~mol.types.mask] = 2
            elif mol.properties['topology'] == 'loop':
              mol.types[~mol.types.mask] = 3
            elif mol.properties['topology'] == 'tie':
              mol.types[~mol.types.mask] = 4


    if accept and not regrow:
      print 'NON REGROW!!'
      self.engine.viz.clear()
      self.engine.viz.draw_system(bonds=True)
      self.engine.viz.show()
      ist()

    self.rosen_subchain1.reset()
    self.rosen_subchain2.reset()
    self.rosen_chain.reset()
    mc_move_data['string'] += ' Wnew/Wold: {:3.2e}/{:3.2e}={:5.4f}'.format(Wnew,Wold,Wnew/Wold)
    return accept,mc_move_data
  def pick_mol(self):
    tail_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tail']
    tie_list  = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='tie']
    loop_list = [mol for mol in self.system.molecule_types['ChainSegment'] if mol.properties['topology']=='loop']

    if len(tail_list)==0:
      abort = True
      counts = {'tail':0,'tie':0,'loop':0)
      return abort,None,None,counts

    index_list_tail = [i for mol in tail_list for i in mol.indices]
    index_list_tie  = [i for mol in tail_list for i in mol.indices]
    index_list_loop = [i for mol in tail_list for i in mol.indices]
    counts = {'tail':len(index_list_tail),'tie':len(index_list_tie),'loop':len(index_list_loop))

    index_list, = np.where(np.in1d(self.system.types,self.regrowth_types))
    index1 = choice(index_list)
    index_list[index1] = False
    mol1 = self.system.molecule_map[index1]

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
  def set_indices(self,chain1_dict,chain2_dict):
    indices1 = chain1_dict['mol'].indices
    indices2 = chain2_dict['mol'].indices

    # We want to be growing in the direction away from the crystallite
    if indices1[-1] in chain1_dict['mol'].properties['connected_to']:
      indices1 = indices1[::-1]

    if indices2[-1] in chain2_dict['mol'].properties['connected_to']:
      indices2 = indices2[::-1]

    self.rosen_subchain1.set_indices(indices1,internal_growth=False,random_inversion=False)
    self.rosen_subchain2.set_indices(indices2,internal_growth=False,random_inversion=False)

    self.rosen_subchain1.build_arrays()
    self.rosen_subchain2.build_arrays(index_shift  = self.rosen_subchain1.length)

    # The end of chain2 will be the beacon of chain1 and vice versa
    chain1_end_sys = self.rosen_subchain1.indices[-1]
    chain2_end_sys = self.rosen_subchain2.indices[-1]
    chain1_start_sys = self.rosen_subchain1.indices[0]
    chain2_start_sys = self.rosen_subchain2.indices[0]
    chain2_start_outer_sys = self.rosen_subchain2.retrace_bonds[0][0][0]
    chain1_start_outer_sys = self.rosen_subchain1.retrace_bonds[0][0][0]

    l1 = self.rosen_subchain1.length
    l2 = self.rosen_subchain2.length
    self.rosen_subchain1.regrowth_beacons = [(chain2_start_outer_sys,l1+l2-i) for i in range(l1)] 
    self.rosen_subchain2.regrowth_beacons = [(chain2_start_outer_sys,l2-i) for i in range(l2)]
    self.rosen_subchain1.regrowth_anchors = self.rosen_subchain1.retrace_anchors[:]
    self.rosen_subchain2.regrowth_anchors = [None for _ in self.rosen_subchain2.indices]

    # old chain: 0--1--2--3--4--5--6--7--8--9--10
    # new chain: 0--1--2--3--4--5  6--7--8--9--10
    # growing the old chain (steps)
    #      x-------------->
    #   1) 0--1--2--3--4--5
    #                      <--------------x
    #   2) 0--1--2--3--4--5--6--7--8--9--10
    #   i.e. 5-6 bond gets added from "chain2"
    # ...this is a PITA and is probably wrong in some way
    self.rosen_subchain1.regrowth_bonds = copy.deepcopy(self.rosen_subchain1.retrace_bonds)
    self.rosen_subchain2.regrowth_bonds = copy.deepcopy(self.rosen_subchain2.retrace_bonds)


    chain1_end_trial = self.rosen_subchain1.regrowth_bonds[-1][-1][-1]
    chain2_end_trial = self.rosen_subchain2.regrowth_bonds[-1][-1][-1]
    self.rosen_subchain2.regrowth_bonds[-1].append([chain2_end_trial,chain2_start_outer_sys])
    self.rosen_subchain2.regrowth_bonds[0][-1][0] = chain1_end_trial

    ## Combine the two rosenbluth chains into one
    self.rosen_chain.combine(self.rosen_subchain1,self.rosen_subchain2)

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([chain1_end_sys,chain2_start_sys])
    bonds['added'].append([chain2_end_sys,chain2_start_outer_sys])
    bonds['removed'].append([chain2_start_sys,chain2_start_outer_sys])

    mol1 = {}
    mol1['molecule']     = chain1_dict['mol']
    mol1['indices'] =  self.rosen_subchain1.low_indices + self.rosen_subchain1.indices
    mol1['indices'] += self.rosen_subchain2.indices + self.rosen_subchain2.low_indices[::-1]

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'mol1':mol1,'added':[],'removed':[chain2_dict['mol']]}

    # print "CHAIN1:"
    # print chain1_dict['mol'].indices
    # print "CHAIN2:"
    # print chain2_dict['mol'].indices
    # print "GROW1:"
    # print self.rosen_subchain1.indices
    # print "GROW2:"
    # print self.rosen_subchain2.indices
    # print "FULL BONDS"
    # print self.rosen_chain.regrowth_bonds
    # print "BOND DELTA"
    # print bonds
    # ist() # splice
  def set_indices1(self,chain1_dict,chain2_dict):
    indices1 = chain1_dict['mol'].indices
    indices2 = chain2_dict['mol'].indices

    # We want to be growing in the direction away from the crystallite
    if indices1[-1] in chain1_dict['mol'].properties['connected_to']:
      indices1 = indices1[::-1]

    if indices2[-1] not in chain2_dict['mol'].properties['connected_to']:
      indices2 = indices2[::-1]

    self.rosen_subchain1.set_indices(indices1,internal_growth=False,random_inversion=False)
    self.rosen_subchain2.set_indices(indices2,internal_growth=False,random_inversion=False)

    self.rosen_subchain1.build_arrays()
    self.rosen_subchain2.build_arrays(index_shift  = self.rosen_subchain1.length)

    # The end of chain2 will be the beacon of chain1 and vice versa
    chain1_end_sys = self.rosen_subchain1.indices[-1]
    chain2_end_sys = self.rosen_subchain2.indices[-1]
    chain1_start_sys = self.rosen_subchain1.indices[0]
    chain2_start_sys = self.rosen_subchain2.indices[0]
    chain2_start_outer_sys = self.rosen_subchain2.retrace_bonds[0][0][0]
    chain1_start_outer_sys = self.rosen_subchain1.retrace_bonds[0][0][0]

    l1 = self.rosen_subchain1.length
    l2 = self.rosen_subchain2.length
    self.rosen_subchain1.regrowth_beacons = [(indices2[0],l1-i) for i in range(l1)] 
    self.rosen_subchain2.regrowth_beacons = self.rosen_subchain2.retrace_beacons[:]
    self.rosen_subchain1.regrowth_anchors = self.rosen_subchain1.retrace_anchors[:]
    self.rosen_subchain2.regrowth_anchors = self.rosen_subchain2.retrace_anchors[:]

    # old chain: 0--1--2--3--4--5--6--7--8--9--10
    # new chain: 0--1--2--3--4--5  6--7--8--9--10
    # growing the old chain (steps)
    #      x-------------->
    #   1) 0--1--2--3--4--5
    #                      <--------------x
    #   2) 0--1--2--3--4--5--6--7--8--9--10
    #   i.e. 5-6 bond gets added from "chain2"
    # ...this is a PITA and is probably wrong in some way
    self.rosen_subchain1.regrowth_bonds = copy.deepcopy(self.rosen_subchain1.retrace_bonds)
    self.rosen_subchain2.regrowth_bonds = copy.deepcopy(self.rosen_subchain2.retrace_bonds)


    chain1_end_trial = self.rosen_subchain1.regrowth_bonds[-1][-1][-1]
    self.rosen_subchain1.regrowth_bonds[-1].append([chain1_end_trial,indices2[0]])

    ## Combine the two rosenbluth chains into one
    self.rosen_chain.combine(self.rosen_subchain1,self.rosen_subchain2)

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([indices1[-1],indices2[0]])

    mol1 = {}
    mol1['molecule']     = chain1_dict['mol']
    mol1['indices'] =  [i for submol in self.rosen_chain.all_indices for i in submol]

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'mol1':mol1,'added':[],'removed':[chain2_dict['mol']]}

    # print "CHAIN1:"
    # print chain1_dict['mol'].indices
    # print "CHAIN2:"
    # print chain2_dict['mol'].indices
    # print "GROW1:"
    # print self.rosen_subchain1.indices
    # print "GROW2:"
    # print self.rosen_subchain2.indices
    # print "FULL BONDS"
    # print self.rosen_chain.regrowth_bonds
    # print "BOND DELTA"
    # print bonds
    # ist() # splice
  def set_indices2(self,chain1_dict,chain2_dict):
    indices1 = chain1_dict['mol'].indices
    indices2 = chain2_dict['mol'].indices

    # We want to be growing in the direction away from the crystallite
    if indices1[-1] in chain1_dict['mol'].properties['connected_to']:
      indices1 = indices1[::-1]

    if indices2[-1] in chain2_dict['mol'].properties['connected_to']:
      indices2 = indices2[::-1]

    self.rosen_subchain1.set_indices(indices1,internal_growth=False,random_inversion=False)
    self.rosen_subchain2.set_indices(indices2,internal_growth=False,random_inversion=False)

    self.rosen_subchain1.build_arrays()
    self.rosen_subchain2.build_arrays(index_shift  = self.rosen_subchain1.length)

    # The end of chain2 will be the beacon of chain1 and vice versa
    chain1_end_sys = self.rosen_subchain1.indices[-1]
    chain2_end_sys = self.rosen_subchain2.indices[-1]
    chain2_start_sys = self.rosen_subchain2.retrace_bonds[0][0][0]
    chain1_end_trial = self.rosen_subchain1.retrace_bonds[-1][-1][-1]

    l1 = self.rosen_subchain1.length
    l2 = self.rosen_subchain2.length
    self.rosen_subchain1.regrowth_beacons = [(chain2_start_sys,l1+l2-i) for i in range(l1)] 
    # self.rosen_subchain1.regrowth_beacons = [None for i in range(l1)] 
    self.rosen_subchain2.regrowth_beacons = [(chain1_end_trial,l2-i) for i in range(l2)]
    self.rosen_subchain1.regrowth_anchors = self.rosen_subchain1.retrace_anchors[:]
    self.rosen_subchain2.regrowth_anchors = self.rosen_subchain2.retrace_anchors[:]

    # old chain: 0--1--2--3--4--5--6--7--8--9--10
    # new chain: 0--1--2--3--4--5  6--7--8--9--10
    # growing the old chain (steps)
    #      x-------------->
    #   1) 0--1--2--3--4--5
    #                      <--------------x
    #   2) 0--1--2--3--4--5--6--7--8--9--10
    #   i.e. 5-6 bond gets added from "chain2"
    # ...this is a PITA and is probably wrong in some way
    self.rosen_subchain1.regrowth_bonds = copy.deepcopy(self.rosen_subchain1.retrace_bonds)
    self.rosen_subchain2.regrowth_bonds = copy.deepcopy(self.rosen_subchain2.retrace_bonds)


    # Add in the connecting bond
    bondi = self.rosen_subchain1.regrowth_bonds[-1][-1][-1]
    bondj = self.rosen_subchain2.regrowth_bonds[-1][-1][-1]
    self.rosen_subchain2.regrowth_bonds[-1].append([bondi,bondj])

    ## Combine the two rosenbluth chains into one
    self.rosen_chain.combine(self.rosen_subchain1,self.rosen_subchain2)

    ## Need to build `acceptance package` i.e. what needs to be changed (beyond x,y,z,imx,imy,imz) if
    ## the move is accepted. 
    bonds = {'added':[], 'removed':[]}
    bonds['added'].append([chain1_end_sys,chain2_end_sys])

    mol1 = {}
    mol1['molecule']     = chain1_dict['mol']
    mol1['indices'] = self.rosen_chain.indices

    self.rosen_chain.acceptance_package['bonds'] = bonds
    self.rosen_chain.acceptance_package['molecules'] = {'mol1':mol1,'added':[],'removed':[chain2_dict['mol']]}

    # print "RETRACE_BONDS:"
    # print self.rosen_subchain1.retrace_bonds
    # print "REGROWTH_BONDS:"
    # print self.rosen_subchain1.regrowth_bonds
    # ist() # splice
