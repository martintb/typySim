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
  XXX_local_index
  XXX_sys_index
  '''
  def __init__(self,regrowth_types,num_trials=6,regrowth_min=-1, regrowth_max=-1):
    super(CBMCMove,self).__init__() #must call parent class' constructor
    self.name='CBMCMove'
    self.regrowth_types = regrowth_types
    self.num_trials= num_trials
    self.trial_indices= np.arange(num_trials,dtype=np.int)
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max
  @MonteCarloMove.counter
  def attempt(self):
    UOld = self.engine.TPE_list[-1]

    # Choose a chain molecule and a starting point
    mol = self.system.select.random_molecule(names=['ChainSegment'])
    indices = sorted(list(mol.indices))
    chain_length = len(indices)

    if random()>0.5:
      growing_up = True
    else:
      growing_up = False

    if self.regrowth_min == self.regrowth_max == -1:
      if growing_up:
        start_index_local = 0
        end_index_local = chain_length-1
      else:
        start_index_local = chain_length-1
        end_index_local = 0
      start_index_sys = indices[start_index_local]
      end_index_sys = indices[end_index_local]
      regrowth_indices = indices
      regrowth_size = chain_length
    else:
      start_index_local = choice(range(chain_length))
      start_index_sys = indices[start_index_local]

      # Decide is we are growing with increasing indices or decreasing indices
      # XXX The reliance on continuous indices is worrisome. An alternative is needed. 
      if growing_up:
        segment_size = len(indices[start_index_local:])
      else:
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
      raise ValueError('This move only supports linear polymers')

    URegrow = sum(self.engine.TPE.compute(partial_indices=regrowth_indices))
    UBase = UOld-URegrow


    rosen_weights_new,trial_data = self.rosenbluth(UBase,start_index_sys,start_bond,end_bond,regrowth_indices,retrace=False)

    if trial_data['immediate_reject']:
      return False,-2

    rosen_weights_old,_          = self.rosenbluth(UBase,start_index_sys,start_bond,end_bond,regrowth_indices,retrace=True)

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

    if accept:
      for local_index,sys_index in enumerate(regrowth_indices):
        self.system.x[sys_index] = trial_data['x'][local_index]
        self.system.y[sys_index] = trial_data['y'][local_index]
        self.system.z[sys_index] = trial_data['z'][local_index]
      U = sum(self.engine.TPE.compute())
    else:
      U = -3
    return accept,U
  def rosenbluth(self,UBase,start_index_sys,start_bond,end_bond,regrowth_indices,retrace=False):
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
    trial_x = np.empty((num_trials,regrowth_size),dtype=np.float)
    trial_y = np.empty((num_trials,regrowth_size),dtype=np.float)
    trial_z = np.empty((num_trials,regrowth_size),dtype=np.float)
    trial_types = np.empty((num_trials,regrowth_size),dtype=np.int)
    trial_bonds = np.array(start_bond,dtype=np.int)

    rosen_weights = []
    for local_index,sys_index in enumerate(regrowth_indices):
      pertubations = np.random.random((num_trials,3)) - 0.5
      pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T
      trial_potential_energies = []
      for trial_index,vec in enumerate(pertubations):
        new_x = prev_x + vec[0]
        new_y = prev_y + vec[1]
        new_z = prev_z + vec[2]
        new_x, new_y, new_z = self.system.box.numpy_wrap_position(x=new_x, y=new_y, z=new_z)
        trial_x[trial_index,local_index] = new_x
        trial_y[trial_index,local_index] = new_y
        trial_z[trial_index,local_index] = new_z
        trial_types[trial_index,local_index] = self.system.types[sys_index]

      self.system.set_trial_move(
          x=trial_x[:,:(local_index+1)],
          y=trial_y[:,:(local_index+1)],
          z=trial_z[:,:(local_index+1)],
          types=trial_types[:,:(local_index+1)],
          bonds=trial_bonds,
          )
      #calculation will ignore regrowth_indices
      UList = self.engine.TPE.compute( trial_move=True, partial_indices=regrowth_indices, ntrials=num_trials)
      UList = np.sum(np.add(UList,UBase),axis=0)
      trial_potential_energies.extend(UList)

      if retrace:
        trial_x[0,local_index] = self.system.x[sys_index]
        trial_y[0,local_index] = self.system.y[sys_index]
        trial_z[0,local_index] = self.system.z[sys_index]
        trial_types[0,local_index] = self.system.types[sys_index]
        self.system.set_trial_move(
            x=trial_x[:,:(local_index+1)],
            y=trial_y[:,:(local_index+1)],
            z=trial_z[:,:(local_index+1)],
            types=trial_types[:,:(local_index+1)],
            bonds=trial_bonds
            )
        UList = self.engine.TPE.compute(trial_move=True,partial_indices=regrowth_indices,ntrials=1)
        UList = np.sum(np.add(UList,UBase),axis=0)
        trial_potential_energies.append(UList)
        rosen_weights.append(np.exp(np.negative(trial_potential_energies)))
        chosen_index = 0
      else:
        rosen_weights.append(np.exp(np.negative(trial_potential_energies)))
        # print 'i,ROSEN',local_index,rosen_weights[-1]
        if np.sum(rosen_weights[-1])<1e-16:
          trial_data['immediate_reject'] = True
          break
        else:
          trial_probabilities = rosen_weights[-1]/np.sum(rosen_weights[-1])
          chosen_index = choice(self.trial_indices,p=trial_probabilities)

      # set chosen x,y,z,types 
      for trial_index,vec in enumerate(pertubations):
        trial_x[trial_index,local_index] = trial_x[chosen_index,local_index]
        trial_y[trial_index,local_index] = trial_y[chosen_index,local_index]
        trial_z[trial_index,local_index] = trial_z[chosen_index,local_index]


      if len(trial_bonds)==0:
        trial_bonds = np.array([[trial_sys_index+local_index,trial_sys_index+local_index+1]],dtype=np.int)
      else:
        trial_bonds = np.append(trial_bonds,
                                [[trial_sys_index+local_index,trial_sys_index+local_index+1]]
                                ,axis=0).astype(np.int)

      if local_index==(regrowth_size-2) and len(end_bond)>0:
        trial_bonds = np.append(trial_bonds,end_bond,axis=0).astype(np.int)

      prev_x = trial_x[chosen_index,local_index]
      prev_y = trial_y[chosen_index,local_index]
      prev_z = trial_z[chosen_index,local_index]
    trial_data['x'] = trial_x[0]
    trial_data['y'] = trial_y[0]
    trial_data['z'] = trial_z[0]
    return rosen_weights,trial_data









