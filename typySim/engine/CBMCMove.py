from MonteCarloMove import *
from typySim import molecule
from typySim.geometry import linalg
import numpy as np
from numpy.random import randint
from numpy.random import random,choice

import ipdb as pdb; 
import os; ex = lambda : os._exit(1)
ist = pdb.set_trace


class CBMCMove(MonteCarloMove):
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
  def __init__(self,regrowth_types,num_trials=10):
    super(CBMCMove,self).__init__() #must call parent class' constructor
    self.name='CBMCMove'
    self.regrowth_types = regrowth_types
    self.num_trials= num_trials
    self.trial_indices= np.arange(num_trials,dtype=np.int)
    self.regrowth_indices = None
  @MonteCarloMove.counter
  def attempt(self):
    mc_move_data = {}

    ################################
    ## DETERMINE REGROWTH INDICES##
    ################################
    # Choose a chain molecule and a starting point
    regrowth_mol = self.system.select.random_molecule(names=['ChainSegment'])
    indices = sorted(list(regrowth_mol.indices))
    chain_length = len(indices)

    if chain_length<3:
      accept = False
      mc_move_data['string'] = 'chain_is_too_short'
      return accept,mc_move_data

    # We can regrow the chain from the top-down or bottom-up
    if random()>0.5:
      growing_up = True
      # start_index_local = 1
      start_index_local = randint(1,chain_length-3)
      end_index_local = chain_length-1
      self.regrowth_indices  = indices[start_index_local:]
    else:
      growing_up = False
      # start_index_local = chain_length-2
      start_index_local = randint(3,chain_length-1) #start_index_local is [low,high)
      end_index_local = 0
      self.regrowth_indices  = indices[start_index_local::-1]
    start_index_sys        = indices[start_index_local]
    end_index_sys          = indices[end_index_local]
    regrowth_size          = len(self.regrowth_indices)
    start_index_new        = self.system.nbeads
    end_index_new          = self.system.nbeads + abs(regrowth_size)-1

    # Need to analyze bonds to get starting and ending bonds
    start_bond = self.get_outer_bonds(start_index_sys,start_index_new)
    end_bond = self.get_outer_bonds(end_index_sys,end_index_new)

    if len(start_bond)>1 and len(end_bond)>1:
      raise ValueError('This move only supports linear polymers')

    ################################
    ## INITIAL ENERGY CALCULATION ##
    ################################
    # Need to calculate the base energy of the system that we are regrowing into
    UOld    = self.engine.TPE_list[-1]
    URegrow = self.engine.TPE.compute(partial_indices=self.regrowth_indices)
    UBase   = UOld-sum(URegrow)

    ####################
    ## GROW NEW CHAIN ##
    ####################
    rosen_weights_new,trial_data = self.rosenbluth(UBase,start_index_sys,start_bond,end_bond,retrace=False)
    if trial_data['abort']:
      mc_move_data['string'] = 'bad_trial_move_new'
      accept = False
      return accept,mc_move_data

    #####################
    ## TRACE OLD CHAIN ##
    #####################
    rosen_weights_old,old_trial_data = self.rosenbluth(UBase,start_index_sys,start_bond,end_bond,retrace=True)
    if old_trial_data['abort']:
      mc_move_data['string'] = 'bad_trial_move_old'
      accept = False
      return accept,mc_move_data

    #################################
    ## CALCULATE FULL ROSEN FACTOR ##
    #################################
    Wnew = np.product(np.sum(rosen_weights_new,axis=1))
    Wold = np.product(np.sum(rosen_weights_old,axis=1))
    mc_move_data['Wnew'] = Wnew
    mc_move_data['Wold'] = Wold
    mc_move_data['k'] = regrowth_size

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
      for local_index,sys_index in enumerate(self.regrowth_indices):
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
  def get_outer_bonds(self,sys_index,new_index):
    bonds = []
    for j in np.array(self.system.bonds.bonds[sys_index]):
      if j==-1:
        break
      elif j not in self.regrowth_indices:
        bonds.append([j,new_index])
    return bonds
  def rosenbluth(self,UBase,start_index_sys,start_bond,end_bond,retrace=False):
    trial_data = {}
    trial_data['abort'] = False

    regrowth_size = len(self.regrowth_indices)
    start_index_new = self.system.nbeads

    # XXX This isn't a good choice for the case of regrowing the end of the chain
    #     Need to come up with something better
    if len(start_bond)==1:
      prev_index_sys = start_bond[0][0]
    else:
      prev_index_sys = start_index_sys


    prev_x = self.system.x[prev_index_sys]
    prev_y = self.system.y[prev_index_sys]
    prev_z = self.system.z[prev_index_sys]
    trial_x = np.empty((self.num_trials,regrowth_size),dtype=np.float)
    trial_y = np.empty((self.num_trials,regrowth_size),dtype=np.float)
    trial_z = np.empty((self.num_trials,regrowth_size),dtype=np.float)

    prev_imx = self.system.imx[prev_index_sys]
    prev_imy = self.system.imy[prev_index_sys]
    prev_imz = self.system.imz[prev_index_sys]
    trial_imx = np.empty((self.num_trials,regrowth_size),dtype=np.int)
    trial_imy = np.empty((self.num_trials,regrowth_size),dtype=np.int)
    trial_imz = np.empty((self.num_trials,regrowth_size),dtype=np.int)

    trial_types = np.empty((self.num_trials,regrowth_size),dtype=np.int)
    trial_bonds = np.array(start_bond,dtype=np.int)

    rosen_weights = []
    for local_index,sys_index in enumerate(self.regrowth_indices):

      pertubations = np.random.random((self.num_trials,3)) - 0.5
      pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T

      # If this is a retracing move, set the first trial position to be the 
      # actual position of the bead and then remove one of the trial pertubations.
      if retrace:
        trial_x[0,local_index] = self.system.x[sys_index]
        trial_y[0,local_index] = self.system.y[sys_index]
        trial_z[0,local_index] = self.system.z[sys_index]
        trial_imx[0,local_index] = self.system.imx[sys_index]
        trial_imy[0,local_index] = self.system.imy[sys_index]
        trial_imz[0,local_index] = self.system.imz[sys_index]
        trial_types[0,local_index] = self.system.types[sys_index]
        pertubations = pertubations[1:]
        trial_index_start = 1
      else:
        trial_index_start = 0

      for trial_index,vec in enumerate(pertubations,start=trial_index_start):
        new_pos = vec[np.newaxis].T
        new_pos[0] += prev_x
        new_pos[1] += prev_y
        new_pos[2] += prev_z
        (new_x, new_y, new_z),(new_imx,new_imy,new_imz) = self.system.box.wrap_positions(new_pos[0],new_pos[1],new_pos[2])
        trial_x[trial_index,local_index] = new_x[0]
        trial_y[trial_index,local_index] = new_y[0]
        trial_z[trial_index,local_index] = new_z[0]
        trial_imx[trial_index,local_index] = (prev_imx + new_imx[0])
        trial_imy[trial_index,local_index] = (prev_imy + new_imy[0])
        trial_imz[trial_index,local_index] = (prev_imz + new_imz[0])
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

      # If all of the rosen_weights are extremely small or zero, it means that
      # the configuration is "stuck" and that all trial monomers are high energy.
      # We abort the growth early in order to not waste time continuing the growth
      # of a broken configuration.
      # if (not retrace) and np.sum(rosen_weights[-1])<1e-16:
      if (not retrace) and np.all(rosen_weights[-1]==0.):
        trial_data['abort'] = True
        return rosen_weights,trial_data

      if retrace:
        chosen_index = 0
      else:
        trial_probabilities = rosen_weights[-1]/np.sum(rosen_weights[-1])
        chosen_index = choice(self.trial_indices,p=trial_probabilities)

      # set chosen x,y,z,types 
      for trial_index,vec in enumerate(pertubations):
        trial_x[trial_index,local_index] = trial_x[chosen_index,local_index]
        trial_y[trial_index,local_index] = trial_y[chosen_index,local_index]
        trial_z[trial_index,local_index] = trial_z[chosen_index,local_index]
        trial_imx[trial_index,local_index] = trial_imx[chosen_index,local_index]
        trial_imy[trial_index,local_index] = trial_imy[chosen_index,local_index]
        trial_imz[trial_index,local_index] = trial_imz[chosen_index,local_index]
      prev_x = trial_x[chosen_index,local_index]
      prev_y = trial_y[chosen_index,local_index]
      prev_z = trial_z[chosen_index,local_index]
      prev_imx = trial_imx[chosen_index,local_index]
      prev_imy = trial_imy[chosen_index,local_index]
      prev_imz = trial_imz[chosen_index,local_index]

      if len(trial_bonds)==0:
        trial_bonds = np.array([[start_index_new+local_index,start_index_new+local_index+1]],dtype=np.int)
      else:
        trial_bonds = np.append(trial_bonds,
                                [[start_index_new+local_index,start_index_new+local_index+1]]
                                ,axis=0).astype(np.int)

      if local_index==(regrowth_size-2) and len(end_bond)>0:
        trial_bonds = np.append(trial_bonds,end_bond,axis=0).astype(np.int)


    # It doesn't matter which trial_index is chosen as all are equalized at this point. 
    trial_data['x'] = trial_x[chosen_index,:] 
    trial_data['y'] = trial_y[chosen_index,:]
    trial_data['z'] = trial_z[chosen_index,:]
    trial_data['imx'] = trial_imx[chosen_index,:] 
    trial_data['imy'] = trial_imy[chosen_index,:]
    trial_data['imz'] = trial_imz[chosen_index,:]
    trial_data['U'] = trial_potential_energies[chosen_index]

    return rosen_weights,trial_data









