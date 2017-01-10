from MonteCarloMove import *
from typySim import molecule
from typySim.geometry import linalg
import numpy as np
from numpy.random import randint
from numpy.random import random,choice

import ipdb as pdb; 
import os; ex = lambda : os._exit(1)
ist = pdb.set_trace

class InteriorCBMCMove(MonteCarloMove):
  '''
  Definitions
  -----------
  XXX_local_index
  XXX_sys_index
  '''
  def __init__(self,regrowth_types,num_trials=-1,regrowth_min = -1, regrowth_max=6):
    super(InteriorCBMCMove,self).__init__() #must call parent class' constructor
    self.name='InteriorCBMCMove'
    self.regrowth_types = regrowth_types
    self.num_trials= num_trials
    self.trial_indices= np.arange(num_trials,dtype=np.int)
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max
  @MonteCarloMove.counter
  def attempt_old(self):
    Uold = self.engine.TPE_list[-1]
    print 'UOLD',Uold,'NBEADS',self.system.nbeads

    # Choose a chain molecule and a starting point
    mol = self.system.select.random_molecule(names=['ChainSegment'])
    indices = sorted(list(mol.indices))
    chain_length = len(indices)
    start_index_local = choice(range(chain_length))
    start_index_sys = indices[start_index_local]
    print 'indices',indices

    # Decide is we are growing with increasing indices or decreasing indices
    # XXX The reliance on continuous indices is worrisome. An alternative is needed. 
    if random()>0.5:
      growing_up = True
      segment_size = len(indices[start_index_local:])
    else:
      growing_up = False
      segment_size = len(indices[:start_index_local])

    # If there aren't enough beads to regrow, reject the move
    if segment_size<self.regrowth_min:
      return False,-1

    print 'N',chain_length,'start_index_local/sys',start_index_local,start_index_sys,'growing_up=',growing_up

    # Determing the starting and ending index of regrowth
    if growing_up:
      minsize = self.regrowth_min
      maxsize = min([segment_size,self.regrowth_max])
      if maxsize==minsize:
        regrowth_size = maxsize
      else:
        regrowth_size = randint(minsize,maxsize)
      end_index_local = start_index_local + regrowth_size-1
      end_index_sys = indices[end_index_local]
      regrowth_indices = indices[start_index_local:(end_index_local+1)]
    else:
      minsize = -self.regrowth_min
      maxsize = -min([segment_size,self.regrowth_max])
      if maxsize==minsize:
        regrowth_size = maxsize
      else:
        regrowth_size = randint(maxsize,minsize)
      end_index_local = start_index_local + regrowth_size+1
      end_index_sys = indices[end_index_local]
      regrowth_indices = indices[start_index_local:(end_index_local-1):-1]
    trial_sys_index = self.system.nbeads


    print 'min/maxval',minsize,maxsize,'end_index_local/sys',end_index_local,end_index_sys
    print 'Regrowth_indices',regrowth_indices,regrowth_size

    # Need to analyze bonds to get starting and ending bonds
    start_bond = []
    start_bondlist = np.array(self.system.bonds.bonds[start_index_sys])
    for idex in start_bondlist:
      if idex==-1:
        break
      elif idex not in regrowth_indices:
        start_bond.append([idex,trial_sys_index]) #the trial move index is not the same as the start_sys_index

    end_bond = []
    end_bondlist = np.array(self.system.bonds.bonds[end_index_sys])
    for idex in end_bondlist:
      if idex==-1:
        break
      elif idex not in regrowth_indices:
        end_bond.append([trial_sys_index+abs(regrowth_size),idex]) 

    if len(start_bond)>1 and len(end_bond)>1:
      ist()
      raise ValueError('This move only supports linear polymers')

    # Uregrow = self.engine.TPE.compute(partial_indices=regrowth_indices)
    print 'START_BOND',start_bondlist,start_bond
    print 'END_BOND',end_bondlist,end_bond


    print '================ NEW GROW ================'
    ######################
    ## GROW NEW SECTION ##
    ######################
    if len(start_bond)==1:
      prev_index_sys = start_bond[0][0]
    else:
      prev_index_sys = start_index_sys
    prev_x = self.system.x[prev_index_sys]
    prev_y = self.system.y[prev_index_sys]
    prev_z = self.system.z[prev_index_sys]
    new_x = np.array([0],dtype=np.float)
    new_y = np.array([0],dtype=np.float)
    new_z = np.array([0],dtype=np.float)
    new_types = np.array([0],dtype=np.int)
    new_bonds = np.array(start_bond,dtype=np.int)
    rosen_weights_new = []
    for local_index,sys_index in enumerate(regrowth_indices):
      print new_bonds
      pertubations = np.random.random((self.num_trials,3)) - 0.5
      pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T
      Ulist = []
      for vec in pertubations:
        new_x[local_index] = prev_x + vec[0]
        new_y[local_index] = prev_y + vec[1]
        new_z[local_index] = prev_z + vec[2]
        new_types[local_index] = self.system.types[sys_index]
        self.system.set_trial_move(
                                   x=new_x,
                                   y=new_y,
                                   z=new_z,
                                   types=new_types,
                                   bonds=new_bonds
                                  )
        #calculation will ignore regrowth_indices
        Ulist.append(self.engine.TPE.compute(trial_move=True,partial_indices=regrowth_indices))
      # rosen_weights = np.exp(np.negative(np.add(Ulist,Uold-Uregrow)))
      trial_rosen_weights = np.exp(np.negative(Ulist))
      trial_probabilities = trial_rosen_weights/np.sum(trial_rosen_weights)
      chosen_index = choice(self.trial_indices,p=trial_probabilities)

      #store chosen rosen_weight
      rosen_weights_new.append(trial_rosen_weights)

      # set chosen x,y,z,types and expand trial arrays
      new_x[local_index] = prev_x + pertubations[chosen_index,0]
      new_y[local_index] = prev_y + pertubations[chosen_index,1]
      new_z[local_index] = prev_z + pertubations[chosen_index,2]
      new_x = np.append(new_x,[0])
      new_y = np.append(new_y,[0])
      new_z = np.append(new_z,[0])
      new_types = np.append(new_types,[0])

      # update trial_bondlist
      if local_index==(len(regrowth_indices)-1):
        new_bonds = np.append(new_bonds,end_bond,axis=0).astype(np.int)
      else:
        new_bonds = np.append(new_bonds,[[trial_sys_index+local_index,trial_sys_index+local_index+1]],axis=0).astype(np.int)


      prev_x = new_x[local_index]
      prev_y = new_y[local_index]
      prev_z = new_z[local_index]

    print '================ OLD TRACE ================'
    #######################
    ## TRACE OLD SECTION ##
    #######################
    if len(start_bond)==1:
      prev_index_sys = start_bond[0][0]
    else:
      prev_index_sys = start_index_sys
    prev_x = self.system.x[prev_index_sys]
    prev_y = self.system.y[prev_index_sys]
    prev_z = self.system.z[prev_index_sys]
    old_x = np.array([0],dtype=np.float)
    old_y = np.array([0],dtype=np.float)
    old_z = np.array([0],dtype=np.float)
    old_types = np.array([0],dtype=np.int)
    old_bonds = np.array(start_bond,dtype=np.int)
    rosen_weights_old = []
    for local_index,sys_index in enumerate(regrowth_indices):
      print old_bonds
      pertubations = np.random.random((self.num_trials-1,3)) - 0.5
      pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T

      Ulist = []
      for vec in pertubations:
        old_x[local_index] = prev_x + vec[0]
        old_y[local_index] = prev_y + vec[1]
        old_z[local_index] = prev_z + vec[2]
        old_types[local_index] = self.system.types[sys_index]
        self.system.set_trial_move(
                                   x=old_x,
                                   y=old_y,
                                   z=old_z,
                                   types=old_types,
                                   bonds=old_bonds
                                  )
        #calculation will ignore regrowth_indices
        Ulist.append(self.engine.TPE.compute(trial_move=True,partial_indices=regrowth_indices))

      # calculate old PE
      old_x[local_index] = self.system.x[sys_index]
      old_y[local_index] = self.system.y[sys_index]
      old_z[local_index] = self.system.z[sys_index]
      old_types[local_index] = self.system.types[sys_index]
      self.system.set_trial_move(
                                 x=old_x,
                                 y=old_y,
                                 z=old_z,
                                 types=old_types,
                                 bonds=old_bonds
                                )
      Ulist.append(self.engine.TPE.compute(trial_move=True,partial_indices=regrowth_indices))

      trial_rosen_weights = np.exp(np.negative(Ulist))
      rosen_weights_old.append(trial_rosen_weights)

      # set chosen x,y,z,types and expand trial arrays
      old_x = np.append(old_x,[0])
      old_y = np.append(old_y,[0])
      old_z = np.append(old_z,[0])
      old_types = np.append(old_types,[0])

      # update trial_bondlist
      if local_index==(len(regrowth_indices)-1):
        old_bonds = np.append(old_bonds,end_bond,axis=0).astype(np.int)
      else:
        old_bonds = np.append(old_bonds,[[trial_sys_index+local_index,trial_sys_index+local_index+1]],axis=0).astype(np.int)


      prev_x = old_x[local_index]
      prev_y = old_y[local_index]
      prev_z = old_z[local_index]

    Wnew = np.product(np.sum(rosen_weights_new,axis=1))
    Wold = np.product(np.sum(rosen_weights_old,axis=1))
    print 'Wnew/Wold',Wnew,Wold
    ist()
  @MonteCarloMove.counter
  def attempt(self):
    Uold = self.engine.TPE_list[-1]

    # Choose a chain molecule and a starting point
    mol = self.system.select.random_molecule(names=['ChainSegment'])
    indices = sorted(list(mol.indices))
    chain_length = len(indices)

    if self.regrowth_min == self.regrowth_max == -1:
      if random()>0.5:
        start_index_local = 0
        end_index_local = chain_length-1
        growing_up = True
      else:
        start_index_local = chain_length-1
        end_index_local = 0
        growing_up = False
      start_index_sys = indices[start_index_local]
      end_index_sys = indices[end_index_local]
      regrowth_indices = indices
      regrowth_size = chain_length
    else:
      start_index_local = choice(range(chain_length))
      start_index_sys = indices[start_index_local]

      # Decide is we are growing with increasing indices or decreasing indices
      # XXX The reliance on continuous indices is worrisome. An alternative is needed. 
      if random()>0.5:
        growing_up = True
        segment_size = len(indices[start_index_local:])
      else:
        growing_up = False
        segment_size = len(indices[:start_index_local])

      # If there aren't enough beads to regrow, reject the move
      if segment_size<self.regrowth_min:
        return False,-1

      # Determing the starting and ending index of regrowth
      if growing_up:
        minsize = self.regrowth_min
        maxsize = min([segment_size,self.regrowth_max])
        if maxsize==minsize:
          regrowth_size = maxsize
        else:
          regrowth_size = randint(minsize,maxsize)
        end_index_local = start_index_local + regrowth_size-1
        end_index_sys = indices[end_index_local]
        regrowth_indices = indices[start_index_local:(end_index_local+1)]
      else:
        minsize = -self.regrowth_min
        maxsize = -min([segment_size,self.regrowth_max])
        if maxsize==minsize:
          regrowth_size = maxsize
        else:
          regrowth_size = randint(maxsize,minsize)
        end_index_local = start_index_local + regrowth_size+1
        end_index_sys = indices[end_index_local]
        regrowth_indices = indices[start_index_local:(end_index_local-1):-1]
    trial_sys_index = self.system.nbeads


    # Need to analyze bonds to get starting and ending bonds
    start_bond = []
    start_bondlist = np.array(self.system.bonds.bonds[start_index_sys])
    for idex in start_bondlist:
      if idex==-1:
        break
      elif idex not in regrowth_indices:
        start_bond.append([idex,trial_sys_index]) #the trial move index is not the same as the start_sys_index

    end_bond = []
    end_bondlist = np.array(self.system.bonds.bonds[end_index_sys])
    for idex in end_bondlist:
      if idex==-1:
        break
      elif idex not in regrowth_indices:
        end_bond.append([trial_sys_index+abs(regrowth_size)-1,idex]) 

    if len(start_bond)>1 and len(end_bond)>1:
      ist()
      raise ValueError('This move only supports linear polymers')

    # Uregrow = self.engine.TPE.compute(partial_indices=regrowth_indices)
    print 'START_BOND',start_bondlist,start_bond
    print 'END_BOND',end_bondlist,end_bond


    print '================ NEW GROW ================'
    rosen_weights_new,trial_data = self.rosenbluth_sample(start_bond,end_bond,regrowth_indices,retrace=False)

    if trial_data['immediate_reject'] == False:
      return False,-2

    print '================ OLD TRACE ================'
    rosen_weights_old,_          = self.rosenbluth_sample(start_bond,end_bond,regrowth_indices,retrace=True)

    Wnew = np.product(np.sum(rosen_weights_new,axis=1))
    Wold = np.product(np.sum(rosen_weights_old,axis=1))

    if Wnew>Wold:
      accept = True
    else:
      ranf = np.random.rand()
      if ranf<(Wnew/Wold):
        accept=True
      else:
        accept=False

    print 'Wnew,Wold,Wnew/Wold,accept',Wnew,Wold,Wnew/Wold,accept
    ist()
    
    if accept:
      for local_index,sys_index in enumerate(regrowth_indices):
        self.system.x[sys_index] = trial_data['x'][local_index]
        self.system.y[sys_index] = trial_data['y'][local_index]
        self.system.z[sys_index] = trial_data['z'][local_index]
      U = self.engine.TPE.compute()
    else:
      U=-3
    return accept,U

  def rosenbluth_sample(self,start_bond,end_bond,regrowth_indices,retrace=False):
    trial_data = {}
    trial_data['immediate_reject'] = False
    regrowth_size = len(regrowth_indices)

    if retrace:
      num_trials = self.num_trials - 1
    else:
      num_trials = self.num_trials

    # XXX This isn't a good choice for the case of regrowing the end of the chain
    #     Need to come up with something better
    if len(start_bond)==1:
      prev_index_sys = start_bond[0][0]
    else:
      prev_index_sys = start_index_sys

    trial_sys_index = self.system.nbeads

    prev_x = self.system.x[prev_index_sys]
    prev_y = self.system.y[prev_index_sys]
    prev_z = self.system.z[prev_index_sys]
    trial_x = np.empty(regrowth_size,dtype=np.float)
    trial_y = np.empty(regrowth_size,dtype=np.float)
    trial_z = np.empty(regrowth_size,dtype=np.float)
    trial_types = np.empty(regrowth_size,dtype=np.int)
    trial_bonds = np.array(start_bond,dtype=np.int)

    rosen_weights = []
    for local_index,sys_index in enumerate(regrowth_indices):
      pertubations = np.random.random((num_trials,3)) - 0.5
      pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T
      trial_potential_energies = []
      for vec in pertubations:
        trial_x[local_index] = prev_x + vec[0]
        trial_y[local_index] = prev_y + vec[1]
        trial_z[local_index] = prev_z + vec[2]
        trial_types[local_index] = self.system.types[sys_index]
        self.system.set_trial_move(
                                   x=trial_x[:(local_index+1)],
                                   y=trial_y[:(local_index+1)],
                                   z=trial_z[:(local_index+1)],
                                   types=trial_types[:(local_index+1)],
                                   bonds=trial_bonds
                                  )
        #calculation will ignore regrowth_indices
        trial_potential_energies.append(self.engine.TPE.compute(trial_move=True,partial_indices=regrowth_indices))


      if retrace:
        trial_x[local_index] = self.system.x[sys_index]
        trial_y[local_index] = self.system.y[sys_index]
        trial_z[local_index] = self.system.z[sys_index]
        trial_types[local_index] = self.system.types[sys_index]
        self.system.set_trial_move(
                                   x=trial_x[:(local_index+1)],
                                   y=trial_y[:(local_index+1)],
                                   z=trial_z[:(local_index+1)],
                                   types=trial_types[:(local_index+1)],
                                   bonds=trial_bonds
                                  )
        trial_potential_energies.append(self.engine.TPE.compute(trial_move=True,partial_indices=regrowth_indices))
        rosen_weights.append(np.exp(np.negative(trial_potential_energies)))
      else:
        rosen_weights.append(np.exp(np.negative(trial_potential_energies)))
        if np.sum(rosen_weights[-1])<1e-10:
          trial_data['immediate_reject'] = True
          break
        else:
          trial_probabilities = rosen_weights[-1]/np.sum(rosen_weights[-1])
          chosen_index = choice(self.trial_indices,p=trial_probabilities)

          # set chosen x,y,z,types 
          trial_x[local_index] = prev_x + pertubations[chosen_index,0]
          trial_y[local_index] = prev_y + pertubations[chosen_index,1]
          trial_z[local_index] = prev_z + pertubations[chosen_index,2]
      print 'ROSEN',rosen_weights[-1]


      trial_bonds = np.append(trial_bonds,
                              [[trial_sys_index+local_index,trial_sys_index+local_index+1]]
                              ,axis=0).astype(np.int)

      if local_index==(regrowth_size-2) and len(end_bond)>0:
        trial_bonds = np.append(trial_bonds,end_bond,axis=0).astype(np.int)

      prev_x = trial_x[local_index]
      prev_y = trial_y[local_index]
      prev_z = trial_z[local_index]
    trial_data['x'] = trial_x
    trial_data['y'] = trial_y
    trial_data['z'] = trial_z
    return rosen_weights,trial_data









