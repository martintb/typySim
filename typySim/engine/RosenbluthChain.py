from typySim.geometry import linalg
from numpy.random import randint,random,choice
import numpy as np
import ipdb; ist = ipdb.set_trace

import copy

class RosenbluthChain(object):
  ''' Configurational Bias Monte Carlo Helper Class

  This class automagically sets up much of the information a user needs to
  conduct Configurational Bias Monte Carlo simulations (CBMC). It also handles
  the rosenbluth calculation itself.

  .. Note::
      This is not a chain in the "polymer chain" sense, but rather a "chain" of
      beads that are being removed and regrown in a "chain". The beads of this chain
      have no restriction on who or what they are bonded to and could be free beads. 

  Attributes
  ----------
  

  Notes
  -----
  anchor
      The previous position that the current bead is being grown from. This
      could either be the previous bead in the chain, or a point in space. 
      Random pertubations away from this position serve as the possibilites in
      each step of the CBMC process

  beacon
      A position in space that the chain is growing towards. This is only meaningful 
      when conducting Fixed Endpoint CBMC (FE-CBMC). If specified, the user supplied
      guiding bias will be applied based on the distance and number of bonds between
      the current bead and the guide bead.__import__

  regrowth vs retrace
      Regrowth is the process of sequentially 'growing' random new positions for each bead, choosing
      the position based on its Boltzmann weight, and tracking the "partition function" of each step.
      Retrace is the identical process except the "real" position of the bead is included in the list 
      of options and is always chosen. The goal of the retrace step is simply to calculate the "partition
      function" of the original configuration. The comparison of these partition function is what determines
      the acceptance probability of this move. 

  '''
  def __init__(self,num_trials,regrowth_min,regrowth_max,bias,viz=None):
    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max
    self.bias         = bias
    self.reset()
    self.viz = viz
  @property
  def length(self):
    if self.indices is None:
      return None
    else:
      return len(self.indices)
  def reset(self):
    self.internal_growth    = None
    self.growing_up         = None
    self.all_indices        = None

    self.indices            = None
    self.new_indices        = None

    self.retrace_anchors    = None
    self.retrace_beacons    = None
    self.retrace_bonds      = None

    self.regrowth_anchors   = None
    self.regrowth_beacons   = None
    self.regrowth_bonds     = None

    self.chain_type   = ''

    self.acceptance_package = {}
  def combine(self,chain1,chain2):
    self.internal_growth    = [chain1.internal_growth,chain2.internal_growth]
    self.growing_up         = [chain1.growing_up     ,chain2.growing_up]
    self.all_indices        = [chain1.all_indices        ,chain2.all_indices]

    # Needed for rosenbluth calculation therefore flat lists are necessary
    self.new_indices            = chain1.new_indices      + chain2.new_indices
    self.indices                = chain1.indices          + chain2.indices
    self.retrace_anchors        = chain1.retrace_anchors  + chain2.retrace_anchors
    self.retrace_beacons        = chain1.retrace_beacons  + chain2.retrace_beacons
    self.retrace_bonds          = chain1.retrace_bonds    + chain2.retrace_bonds

    if (chain1.regrowth_anchors is not None) and (chain2.regrowth_anchors is not None):
      self.regrowth_anchors = chain1.regrowth_anchors    + chain2.regrowth_anchors
    if (chain1.regrowth_beacons is not None) and (chain2.regrowth_beacons is not None):
      self.regrowth_beacons = chain1.regrowth_beacons    + chain2.regrowth_beacons
    if (chain1.regrowth_bonds is not None) and (chain2.regrowth_bonds is not None):
      self.regrowth_bonds  =  chain1.regrowth_bonds      + chain2.regrowth_bonds
  def copy_retrace_to_regrowth(self):
    self.regrowth_bonds   = copy.deepcopy(self.retrace_bonds)
    self.regrowth_anchors = self.retrace_anchors[:]
    self.regrowth_beacons = self.retrace_beacons[:]
  def set_indices(self,indices,internal_growth,full_chain=False,random_inversion=True):
    self.internal_growth = internal_growth

    if random_inversion and random()>0.5:
      self.growing_up = False
      self.all_indices = indices[::-1]
    else:
      self.growing_up = True
      self.all_indices = indices[::1]

    chain_length = len(self.all_indices)

    if full_chain:
      start_index_local        = 1
      end_index_local          = chain_length-1 
    elif internal_growth:
      start_index_local        = randint(0,chain_length-self.regrowth_min-1)+1
      max_regrowth_index_local = min(start_index_local+self.regrowth_max,chain_length-1) 
      end_index_local     = randint(start_index_local+self.regrowth_min,max_regrowth_index_local+1)-1
    else:
      start_index_local        = randint(0,chain_length-self.regrowth_min-1)+1
      end_index_local = chain_length-1 

    self.indices                = self.all_indices[start_index_local:(end_index_local+1)] #self
    self.low_indices            = self.all_indices[:start_index_local] #self
    self.high_indices           = self.all_indices[(end_index_local+1):] #self
  def get_outer_bonds(self,sys_index,new_index):
    old_bonds = []
    new_bonds = []
    for bond_j in list(self.system.bonds.bonds[sys_index]):
      if bond_j==-1:
        break
      elif bond_j not in self.indices:
        new_bonds.append([bond_j,new_index])
        old_bonds.append([bond_j,sys_index])
    return old_bonds,new_bonds
  def build_arrays(self,index_shift=0):
    last_local_index = (self.length-1)
    ####################
    ## GROWTH_INDICES ##
    ####################
    self.new_indices = []
    for local_index,sys_index in enumerate(self.indices):
      self.new_indices.append(local_index + self.system.nbeads + index_shift)

    ###########
    ## BONDS ##
    ###########
    self.retrace_bonds = [] 
    for local_index,sys_index in enumerate(self.indices):
      new_index = self.new_indices[local_index]
      self.retrace_bonds.append([])

      # Inner Bonds
      if local_index>0:
          self.retrace_bonds[local_index].append([new_index-1,new_index])

      # Outer Bonds
      if (local_index==0) or (local_index==last_local_index):
        outer_bonds_old,outer_bonds_new = self.get_outer_bonds(sys_index,new_index)
        if len(outer_bonds_new)==1:
          self.retrace_bonds[local_index].extend(outer_bonds_new)
        elif len(outer_bonds_new)>1:
          raise ValueError('Not configured for non-linear polymers!')

    #############
    ## ANCHORS ##
    #############
    self.retrace_anchors = [None for _ in self.indices]
    try:
      # [0][0][0] = [first regrowth step][first bond on this step][first bead of first bond]
      self.retrace_anchors[0] = self.retrace_bonds[0][0][0] #should be a sys_index
    except IndexError:
      raise ValueError('No bond at zero index. typySim is not configured to regrow entire chain!')

    #############
    ## BEACONS ##
    #############
    if len(self.retrace_bonds[-1])>1:
      # [-1][1][0] = [last regrowth step][second bond on this step][first bead of first bond]
      # We want the last bond because it should be the external bond
      beacon = self.retrace_bonds[-1][-1][0] #should be a sys_index
      self.retrace_beacons = [(beacon,self.length-i) for i,_ in enumerate(self.indices)]
    else:
      self.retrace_beacons = [None for _ in self.indices]
  def roll_steplist(self,steplist):
    ''' 
    Takes a list that represents sequential additions and rolls it into
    a compelete list of all additions at each step.

    .. Example::

        >>> A = [[1],[2],[3,3.5],[4]]
        >>> B = roll_steplist(A)
        [[1], [1,2], [1,2,3,3.5],[1,2,3,3.5,4]]
    '''
    rolled_list = []
    num_steps = len(steplist)
    for step in range(num_steps):
      chunklist = []
      for chunk in steplist[:(step+1)]:
        chunklist.extend([val for val in chunk])
      rolled_list.append(chunklist)
    return rolled_list
  def draw_trial(self,anchor_x,anchor_y,anchor_z,orig_x,orig_y,orig_z,trial_x,trial_y,trial_z,local_index):
    self.viz.clear()
    x = self.system.x
    y = self.system.y
    z = self.system.z
    t = self.system.types
    self.viz.draw_glyphs(x,y,z,t)

    x=trial_x[0,:(local_index+1)]
    y=trial_y[0,:(local_index+1)]
    z=trial_z[0,:(local_index+1)]
    t=np.ones_like(x)*8
    self.viz.draw_glyphs(x,y,z,t)

    # x=trial_x[0,:local_index]
    # y=trial_y[0,:local_index]
    # z=trial_z[0,:local_index]
    # t=np.ones_like(x)*8
    # self.viz.draw_glyphs(x,y,z,t)

    # x=trial_x[:,local_index]
    # y=trial_y[:,local_index]
    # z=trial_z[:,local_index]
    # t=np.ones_like(x)*9
    # self.viz.draw_glyphs(x,y,z,t)

    x=orig_x
    y=orig_y
    z=orig_z
    t=np.ones_like(x)*10
    self.viz.draw_glyphs(x,y,z,t)

    if anchor_x:
      x=anchor_x
      y=anchor_y
      z=anchor_z
      t=np.ones_like(x)*7
      self.viz.draw_glyphs(x,y,z,t)

    # mols = list(set([self.system.molecule_map[i] for i in self.indices]))
    # com1 = [mol.properties['center_of_mass'] for mol in mols]
    # [mol.reset() for mol in mols]
    # com2 = [mol.properties['center_of_mass'] for mol in mols]
    # print 'COM:',com1,com2
    # if com1:
    #   com = np.array(com1)
    #   x=com[:,0]
    #   y=com[:,1]
    #   z=com[:,2]
    #   t=np.ones_like(x)*6
    #   self.viz.draw_glyphs(x,y,z,t)

    #   com = np.array(com2)
    #   x=com[:,0]
    #   y=com[:,1]
    #   z=com[:,2]
    #   t=np.ones_like(x)*6
    #   self.viz.draw_glyphs(x,y,z,t)

    self.viz.show(blocking=True,resetCamera=True)
  def apply_acceptance_package(self,color_by_topology=True):
    pkg = self.acceptance_package
    #------------------------#
    # UPDATE COORD AND IMAGE #
    #------------------------#
    for local_index,sys_index in enumerate(self.indices):
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

      for mol in pkg['molecules']['added']:
        self.system.add_molecule(mol)
        mol.reset()
        mol.update_properties() #update topology
        if color_by_topology:
          if mol.properties['topology'] == 'tail':
            mol.types[~mol.types.mask] = 2
          elif mol.properties['topology'] == 'loop':
            mol.types[~mol.types.mask] = 3
          elif mol.properties['topology'] == 'tie':
            mol.types[~mol.types.mask] = 4

      for molDict in pkg['molecules']['modded']:
        mol         = molDict['molecule']
        mol.indices = molDict['indices']
        for index in mol.indices:
          self.system.molecule_map[index] = mol

        mol.reset()
        mol.update_properties() #update topology
        if color_by_topology:
          if mol.properties['topology'] == 'tail':
            mol.types[~mol.types.mask] = 2
          elif mol.properties['topology'] == 'loop':
            mol.types[~mol.types.mask] = 3
          elif mol.properties['topology'] == 'tie':
            mol.types[~mol.types.mask] = 4
  def calc_rosebluth(self,UBase,retrace=False):
    abort = False

    trial_indices= np.arange(self.num_trials,dtype=np.int)

    if retrace:
      bonds   = self.roll_steplist(self.retrace_bonds)
      beacons = self.retrace_beacons
      anchors = self.retrace_anchors
    else:
      bonds   = self.roll_steplist(self.regrowth_bonds)
      beacons = self.regrowth_beacons
      anchors = self.regrowth_anchors


    trial_x     = np.empty((self.num_trials,self.length),dtype=np.float)
    trial_y     = np.empty((self.num_trials,self.length),dtype=np.float)
    trial_z     = np.empty((self.num_trials,self.length),dtype=np.float)
    trial_imx   = np.empty((self.num_trials,self.length),dtype=np.int)
    trial_imy   = np.empty((self.num_trials,self.length),dtype=np.int)
    trial_imz   = np.empty((self.num_trials,self.length),dtype=np.int)
    trial_types = np.empty((self.num_trials,self.length),dtype=np.int)

    rosen_weights = []
    bias_weights = []
    for local_index,sys_index in enumerate(self.indices):
      ###########
      ## BONDS ##
      ###########
      trial_bonds = np.array(bonds[local_index],dtype=np.int)

      #############
      ## ANCHORS ##
      #############
      if anchors[local_index] is None:
        # Not specifying the anchor is only valid after the first bead
        if local_index==0:
          raise ValueError('User must supply anchor for first bead!')

        #Use previously accepted growth.
        anchor_x   = trial_x[0,local_index-1] #first index is arbitrary
        anchor_y   = trial_y[0,local_index-1]
        anchor_z   = trial_z[0,local_index-1]
        anchor_imx = trial_imx[0,local_index-1] #first index is arbitrary
        anchor_imy = trial_imy[0,local_index-1]
        anchor_imz = trial_imz[0,local_index-1]
      else:
        anchor_x   = self.system.x[anchors[local_index]]
        anchor_y   = self.system.y[anchors[local_index]]
        anchor_z   = self.system.z[anchors[local_index]]
        anchor_imx = self.system.imx[anchors[local_index]]
        anchor_imy = self.system.imy[anchors[local_index]]
        anchor_imz = self.system.imz[anchors[local_index]]
        
      ##################################
      ## PERTUBATE & CALC PROBABILITY ##
      ##################################
      pertubations = np.random.random((self.num_trials,3)) - 0.5
      pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T

      # If this is a retracing move, set the first trial position to be the 
      # actual position of the bead and then remove one of the trial pertubations.
      if retrace:
        trial_x[0,local_index]     = self.system.x[sys_index]
        trial_y[0,local_index]     = self.system.y[sys_index]
        trial_z[0,local_index]     = self.system.z[sys_index]
        trial_imx[0,local_index]   = self.system.imx[sys_index]
        trial_imy[0,local_index]   = self.system.imy[sys_index]
        trial_imz[0,local_index]   = self.system.imz[sys_index]
        trial_types[0,local_index] = self.system.types[sys_index]
        pertubations = pertubations[1:] #throw away zero pertubation
        trial_index_start = 1
      else:
        trial_index_start = 0

      for trial_index,vec in enumerate(pertubations,start=trial_index_start):
        trial_types[trial_index,local_index] = self.system.types[sys_index]
        new_pos = vec[np.newaxis].T
        new_pos[0] += anchor_x
        new_pos[1] += anchor_y
        new_pos[2] += anchor_z
        (new_x, new_y, new_z),(new_imx,new_imy,new_imz) = self.system.box.wrap_positions(new_pos[0],new_pos[1],new_pos[2])

        trial_x[trial_index,local_index]     = new_x[0]
        trial_y[trial_index,local_index]     = new_y[0]
        trial_z[trial_index,local_index]     = new_z[0]
        trial_imx[trial_index,local_index]   = (anchor_imx + new_imx[0])
        trial_imy[trial_index,local_index]   = (anchor_imy + new_imy[0])
        trial_imz[trial_index,local_index]   = (anchor_imz + new_imz[0])
        trial_types[trial_index,local_index] = self.system.types[sys_index]

      self.system.set_trial_move(
          x=trial_x[:,:(local_index+1)],
          y=trial_y[:,:(local_index+1)],
          z=trial_z[:,:(local_index+1)],
          imx=trial_imx[:,:(local_index+1)],
          imy=trial_imy[:,:(local_index+1)],
          imz=trial_imz[:,:(local_index+1)],
          types=trial_types[:,:(local_index+1)],
          bonds=trial_bonds,
          )

      #calculation will ignore indices and individually calculate the PE for each trial bead
      UList = self.engine.TPE.compute( trial_move=True, partial_indices=self.indices, ntrials=self.num_trials)
      trial_potential_energies = np.add(np.sum(UList,axis=0),UBase)
      rosen_weights.append(np.exp(np.negative(trial_potential_energies)))

      #############
      ## BEACONS ##
      #############
      if beacons[local_index] is not None:
        end_index,nbonds = beacons[local_index]

        if end_index>self.system.nbeads:
          end_index-=self.system.nbeads
          if end_index>=local_index:
            ist()
            raise ValueError('Ill specified beacon. Can\'t use beacon that hasn\'t been grown yet!')
          end_x = trial_x[0,end_index]
          end_y = trial_y[0,end_index]
          end_z = trial_z[0,end_index]
        else:
          end_x = self.system.x[end_index]
          end_y = self.system.y[end_index]
          end_z = self.system.z[end_index]

        ## Calculate guiding bias for all trial positions
        dx = trial_x[:,local_index] - end_x
        dy = trial_y[:,local_index] - end_y
        dz = trial_z[:,local_index] - end_z
        dr = np.array([dx,dy,dz])
        all_dist = self.system.box.wrap_distances(dr[0],dr[1],dr[2])

        biasx = self.bias[nbonds-1]['x']
        biasy = self.bias[nbonds-1]['y']
        biases = np.interp(all_dist,biasx,biasy,left=0,right=0)

        rosen_weights[-1] *= biases
      else:
        biases = np.ones(self.num_trials)

      ############
      ## CHOOSE ##
      ############
      # If all of the rosen_weights are extremely small or zero, it means that
      # the configuration is "stuck" and that all trial monomers are high energy.
      # We abort the growth early in order to not waste time continuing the growth
      # of a broken configuration.
      # if (not retrace) and np.sum(rosen_weights[-1])<1e-16:
      if (not retrace) and np.all(rosen_weights[-1]==0.):
        # if self.viz is not None:
        #   print 'LOCAL/TOTAL:',local_index,'/',self.length-1
        #   print 'CHAIN_TYPE:',self.chain_type
        #   print 'TRIAL PEs:'
        #   print trial_potential_energies
        #   print 'ULIST'
        #   print UList
        #   print 'BIASES'
        #   print biases
        #   orig_x = [self.system.x[i] for i in self.indices]
        #   orig_y = [self.system.y[i] for i in self.indices]
        #   orig_z = [self.system.z[i] for i in self.indices]
        #   anchor_x = [self.system.x[i] for i in anchors if i is not None]
        #   anchor_y = [self.system.y[i] for i in anchors if i is not None]
        #   anchor_z = [self.system.z[i] for i in anchors if i is not None]
        #   self.draw_trial(anchor_x,anchor_y,anchor_z,orig_x,orig_y,orig_z,trial_x,trial_y,trial_z,local_index)
        # ist()
        abort = True
        return abort,None,None

      if retrace:
        chosen_index = 0
      else:
        trial_probabilities = rosen_weights[-1]/np.sum(rosen_weights[-1])
        chosen_index = choice(trial_indices,p=trial_probabilities)

      #need the chosen bias weight for the final acceptance
      bias_weights.append(biases[chosen_index])

      # set chosen x,y,z,types 
      for trial_index,vec in enumerate(pertubations):
        trial_x[trial_index,local_index]   = trial_x[chosen_index,local_index]
        trial_y[trial_index,local_index]   = trial_y[chosen_index,local_index]
        trial_z[trial_index,local_index]   = trial_z[chosen_index,local_index]
        trial_imx[trial_index,local_index] = trial_imx[chosen_index,local_index]
        trial_imy[trial_index,local_index] = trial_imy[chosen_index,local_index]
        trial_imz[trial_index,local_index] = trial_imz[chosen_index,local_index]


    ######################
    ## WRAP IT UP TO GO ##
    ######################
    if not retrace:
      # It doesn't matter which trial_index is chosen as all are equalized at this point. 
      self.acceptance_package['x']   = trial_x[chosen_index,:] 
      self.acceptance_package['y']   = trial_y[chosen_index,:]
      self.acceptance_package['z']   = trial_z[chosen_index,:]
      self.acceptance_package['imx'] = trial_imx[chosen_index,:] 
      self.acceptance_package['imy'] = trial_imy[chosen_index,:]
      self.acceptance_package['imz'] = trial_imz[chosen_index,:]
      self.acceptance_package['U']   = trial_potential_energies[chosen_index]

      # if self.viz is not None and self.chain_type == 'splice':
      #   print 'LOCAL/TOTAL:',local_index,'/',self.length-1
      #   print 'CHAIN_TYPE:',self.chain_type
      #   print 'TRIAL PEs:'
      #   print trial_potential_energies
      #   print 'ULIST'
      #   print UList
      #   print 'BIASES'
      #   print biases
      #   orig_x = [self.system.x[i] for i in self.indices]
      #   orig_y = [self.system.y[i] for i in self.indices]
      #   orig_z = [self.system.z[i] for i in self.indices]
      #   anchor_x = [self.system.x[i] for i in self.retrace_anchors if i is not None]
      #   anchor_y = [self.system.y[i] for i in self.retrace_anchors if i is not None]
      #   anchor_z = [self.system.z[i] for i in self.retrace_anchors if i is not None]
      #   self.draw_trial(anchor_x,anchor_y,anchor_z,orig_x,orig_y,orig_z,trial_x,trial_y,trial_z,local_index)
      #   ist()

    return abort,rosen_weights,bias_weights
