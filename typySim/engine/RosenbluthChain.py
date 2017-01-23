from numpy.random import randint,random,choice
import numpy as np
import ipdb; ist = ipdb.set_trace

class RosenbluthChain(object):
  ''' Configurational Bias Monte Carlo Helper Class

  This class automagically sets up much of the information a user needs to
  conduct Configurational Bias Monte Carlo simulations (CBMC). It also handles
  the rosenbluth calculation itself.

  Definitions
  -----------
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

  '''
  def __init__(self,num_trials,regrowth_min,regrowth_max,bias):
    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max
    self.bias         = bias
    self.trial_indices= np.arange(num_trials,dtype=np.int)
    self.reset()
  @property
  def regrowth_length(self):
    if self.regrowth_indices is None:
      return None
    else:
      return len(self.regrowth_indices)
  @property
  def chain_length(self):
    if self.indices is None:
      return None
    else:
      return len(self.indices)
  def reset(self):
    self.internal_growth    = None
    self.molecule           = None
    self.growing_up         = None
    self.indices            = None

    # needed for rosenbluth calculation
    self.regrowth_indices   = None
    self.old_anchors        = None
    self.old_beacons        = None
    self.old_bonds          = None
  def __iadd__(self,other):
    self.internal_growth    = [self.internal_growth,other.internal_growth]
    self.molecule           = [self.molecule,other.molecule]
    self.growing_up         = [self.growing_up,other.growing_up]
    self.indices            = [self.indices,other.indices]

    # Needed for rosenbluth calculation therefore flat lists are necessary
    self.regrowth_indices   += other.regrowth_indices
    self.old_anchors        += other.old_anchors
    self.old_beacons        += other.old_beacons
    self.old_bonds          += other.old_bonds
  def set_indices(self,mol,internal_growth):
    self.internal_growth = internal_growth
    self.molecule = mol

    if random()>0.5:
      self.growing_up = False
      self.indices = mol.indices[::1]
    else:
      self.growing_up = True
      self.indices = mol.indices[::-1]

    chain_length = len(self.indices)

    #Note that randint is [low,high)
    start_index_local        = randint(0,chain_length-self.regrowth_min-1)+1
    if internal_growth:
      max_regrowth_index_local = min(start_index_local+self.regrowth_max,chain_length-1) 
      end_index_local     = randint(start_index_local+self.regrowth_min,max_regrowth_index_local+1)-1
    else:
      end_index_local = chain_length-1 

    self.regrowth_indices       = self.indices[start_index_local:(end_index_local+1)] #self
  def get_outer_bonds(self,sys_index,new_index):
    bonds = []
    for bond_j in list(self.system.bonds.bonds[sys_index]):
      if bond_j==-1:
        break
      elif bond_j not in self.regrowth_indices:
        bonds.append([bond_j,new_index])
    return bonds
  def build_old_arrays(self,index_shift=0):
    last_local_index = (self.regrowth_length-1)

    ###########
    ## BONDS ##
    ###########
    self.old_bonds          = [] 
    for local_index,sys_index in enumerate(self.regrowth_indices):
      new_index = local_index + self.system.nbeads + index_shift
      if local_index == 0:
        self.old_bonds.append([])
      else:
        self.old_bonds.append(list(self.old_bonds[local_index-1]))

      ## Inner Bonds
      if local_index>0:
          self.old_bonds[local_index].append([new_index-1,new_index])

      ## Outer Bonds
      if (local_index==0) or (local_index==last_local_index):
        outer_bonds = self.get_outer_bonds(sys_index,new_index)
        if len(outer_bonds)==1:
          self.old_bonds[local_index].extend(outer_bonds)
        elif len(outer_bonds)>1:
          raise ValueError('Not configured for non-linear polymers!')

    #############
    ## ANCHORS ##
    #############
    self.old_anchors = [None for _ in self.regrowth_indices]
    try:
      # [0][0][0] = [first regrowth step][first bond on this step][first bead of first bond]
      self.old_anchors[0] = self.old_bonds[0][0][0] #should be a sys_index
    except IndexError:
      raise ValueError('No bond at zero index. typySim is not configured to regrow entire chain!')

    #############
    ## BEACONS ##
    #############
    if self.internal_growth:
      # [last_local_index][1][0] = [last regrowth step][second bond on this step][first bead of first bond]
      # We want the last bond because it should be the external bond
      beacon = self.old_bonds[last_local_index][-1][0] #should be a sys_index
      self.old_beacons = [(beacon,self.regrowth_length-i) for i,_ in enumerate(self.regrowth_indices)]
    else:
      self.old_beacons = [None for _ in self.regrowth_indices]
  def calc_rosenbluth(self,UBase,retrace=False):
    trial_data = {}
    trial_data['abort'] = False

    if retrace:
      bonds   = self.old_bonds
      beacons = self.old_beacons
      anchors = self.old_anchors
    else:
      bonds   = self.new_bonds
      beacons = self.new_beacons
      anchors = self.new_anchors

    trial_x     = np.empty((self.num_trials,self.regrowth_length),dtype=np.float)
    trial_y     = np.empty((self.num_trials,self.regrowth_length),dtype=np.float)
    trial_z     = np.empty((self.num_trials,self.regrowth_length),dtype=np.float)
    trial_imx   = np.empty((self.num_trials,self.regrowth_length),dtype=np.int)
    trial_imy   = np.empty((self.num_trials,self.regrowth_length),dtype=np.int)
    trial_imz   = np.empty((self.num_trials,self.regrowth_length),dtype=np.int)
    trial_types = np.empty((self.num_trials,self.regrowth_length),dtype=np.int)

    rosen_weights = []
    bias_weights = []
    for local_index,sys_index in enumerate(self.regrowth_indices):
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

      #calculation will ignore regrowth_indices and individually calculate the PE for each trial bead
      UList = self.engine.TPE.compute( trial_move=True, partial_indices=self.regrowth_indices, ntrials=self.num_trials)
      trial_potential_energies = np.add(np.sum(UList,axis=0),UBase)
      rosen_weights.append(np.exp(np.negative(trial_potential_energies)))

      #############
      ## BEACONS ##
      #############
      if beacons[local_index] is not None:
        end_index,nbonds = beacons[local_index]

        end_x = self.system.x[end_index]
        end_y = self.system.y[end_index]
        end_z = self.system.z[end_index]

        ## Calculate guiding bias for all trial positions
        dx = trial_x[:,local_index] - end_x
        dy = trial_y[:,local_index] - end_y
        dz = trial_z[:,local_index] - end_z
        dr = np.array([dx,dy,dz])
        all_dist = self.system.box.wrap_distances(dr[0],dr[1],dr[2])

        biasx = self.bias[nbonds]['x']
        biasy = self.bias[nbonds]['y']
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
        trial_data['abort'] = True
        return rosen_weights,trial_data,None

      if retrace:
        chosen_index = 0
      else:
        trial_probabilities = rosen_weights[-1]/np.sum(rosen_weights[-1])
        chosen_index = choice(self.trial_indices,p=trial_probabilities)

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
    # It doesn't matter which trial_index is chosen as all are equalized at this point. 
    trial_data['x'] = trial_x[chosen_index,:] 
    trial_data['y'] = trial_y[chosen_index,:]
    trial_data['z'] = trial_z[chosen_index,:]
    trial_data['imx'] = trial_imx[chosen_index,:] 
    trial_data['imy'] = trial_imy[chosen_index,:]
    trial_data['imz'] = trial_imz[chosen_index,:]
    trial_data['U'] = trial_potential_energies[chosen_index]

    return rosen_weights,trial_data,bias_weights




    
