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
  def __init__(self,num_trials,regrowth_min,regrowth_max,beacon,viz=None):
    self.num_trials   = num_trials
    self.regrowth_min = regrowth_min
    self.regrowth_max = regrowth_max
    self.beacon       = beacon
    self.rotateQ      = linalg.Quaternion()
    self.viz = viz
    self.reset()
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
    self.trial_indices        = None

    self.retrace_anchors    = None
    self.retrace_beacons    = None

    self.regrowth_anchors   = None
    self.regrowth_beacons   = None

    self.chain_type   = ''

    self.acceptance_package = {}
  def combine(self,chain1,chain2):
    self.internal_growth    = [chain1.internal_growth,chain2.internal_growth]
    self.growing_up         = [chain1.growing_up     ,chain2.growing_up]
    self.all_indices        = [chain1.all_indices        ,chain2.all_indices]

    # Needed for rosenbluth calculation therefore flat lists are necessary
    self.trial_indices          = chain1.trial_indices      + chain2.trial_indices
    self.indices                = chain1.indices          + chain2.indices
    self.retrace_anchors        = chain1.retrace_anchors  + chain2.retrace_anchors
    self.retrace_beacons        = chain1.retrace_beacons  + chain2.retrace_beacons

    if (chain1.regrowth_anchors is not None) and (chain2.regrowth_anchors is not None):
      self.regrowth_anchors = chain1.regrowth_anchors    + chain2.regrowth_anchors
    if (chain1.regrowth_beacons is not None) and (chain2.regrowth_beacons is not None):
      self.regrowth_beacons = chain1.regrowth_beacons    + chain2.regrowth_beacons
  def copy_retrace_to_regrowth(self):
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
      end_index_local          = randint(start_index_local+self.regrowth_min,max_regrowth_index_local+1)-1
    else:
      # start_index_local        = randint(0,chain_length-self.regrowth_min-1)+1
      start_index_local        = randint(chain_length-self.regrowth_max,chain_length-self.regrowth_min-1)+1
      end_index_local          = chain_length-1 

    self.indices                = self.all_indices[start_index_local:(end_index_local+1)] #self
    self.low_indices            = self.all_indices[:start_index_local] #self
    self.high_indices           = self.all_indices[(end_index_local+1):] #self
  def set_indices_around(self,indices,around,full_chain=False,random_inversion=True):
    self.internal_growth = True

    if random_inversion and random()>0.5:
      self.growing_up = False
      self.all_indices = indices[::-1]
    else:
      self.growing_up = True
      self.all_indices = indices[::1]

    chain_length = len(self.all_indices)

    start_index_local        = list(self.all_indices).index(around)
    max_regrowth_index_local = min(start_index_local+self.regrowth_max,chain_length-1) 

    bound1 = start_index_local+self.regrowth_min
    bound2 = max_regrowth_index_local+1
    end_index_local       = randint(bound1,bound2)-1

    self.indices                = self.all_indices[start_index_local:(end_index_local+1)] #self
    self.low_indices            = self.all_indices[:start_index_local] #self
    self.high_indices           = self.all_indices[(end_index_local+1):] #self
  def get_outer_anchors(self,sys_index):
    anchor_list = []
    for bond_j in list(self.system.bonds.bonds[sys_index]):
      if bond_j==-1:
        break
      elif bond_j not in self.indices:
        anchor_list.append(bond_j)
    return anchor_list
  def build_arrays(self,index_shift=0,linear_chain=True):
    last_local_index = (self.length-1)

    self.trial_indices = []
    self.retrace_anchors = []
    self.retrace_beacons
    for local_index,sys_index in enumerate(self.indices):
      trial_index = local_index + self.system.nbeads + index_shift
      self.trial_indices.append(trial_index)

      anchors = self.get_outer_anchors(sys_index)
      if len(anchors)>1:
        raise ValueError('RosenbluthChain build error. Too many anchors found.')
      if linear_chain and local_index>0:
        anchors.append(trial_index-1)
      if len(anchors)==0:
        self.retrace_anchors.append(None)
      else:
        self.retrace_anchors.append(tuple(anchors))

    if (anchors is not None) and (len(anchors)>1):
      beacon = anchors[0]
      self.retrace_beacons = [(beacon,self.length-i) for i,_ in enumerate(self.indices)]
    else:
      self.retrace_beacons = [None for _ in self.indices]
  def draw_trial(self,anchors,beacon,grown,trials,orig,chosen):
    self.viz.clear()
    x = self.system.x
    y = self.system.y
    z = self.system.z
    t = self.system.types
    self.viz.draw_glyphs(x,y,z,t)

    x = orig[:,0]
    y = orig[:,1]
    z = orig[:,2]
    t=np.ones_like(x)*9
    self.viz.draw_glyphs(x,y,z,t)

    x = anchors[:,0]
    y = anchors[:,1]
    z = anchors[:,2]
    t=np.ones_like(x)*8
    self.viz.draw_glyphs(x,y,z,t)

    x = trials[:,0]
    y = trials[:,1]
    z = trials[:,2]
    t=np.ones_like(x)*7
    self.viz.draw_glyphs(x,y,z,t)

    if chosen is not None:
      x = [chosen[0]]
      y = [chosen[1]]
      z = [chosen[2]]
      t=np.ones_like(x)*10
      self.viz.draw_glyphs(x,y,z,t)

    if beacon is not None:
      x = [beacon[0]]
      y = [beacon[1]]
      z = [beacon[2]]
      t=np.ones_like(x)*6
      self.viz.draw_glyphs(x,y,z,t)

    if grown is not None:
      x = grown[:,0]
      y = grown[:,1]
      z = grown[:,2]
      t=np.ones_like(x)*5
      self.viz.draw_glyphs(x,y,z,t)

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
            mol.types[~mol.types.mask] = self.system.NonBondedTable.A2N['TAIL']
          elif mol.properties['topology'] == 'loop':
            mol.types[~mol.types.mask] = self.system.NonBondedTable.A2N['LOOP']
          elif mol.properties['topology'] == 'tie':
            mol.types[~mol.types.mask] = self.system.NonBondedTable.A2N['TIE']

      for molDict in pkg['molecules']['modded']:
        mol         = molDict['molecule']
        mol.indices = molDict['indices']
        for index in mol.indices:
          self.system.molecule_map[index] = mol

        mol.reset()
        mol.update_properties() #update topology
        if color_by_topology:
          if mol.properties['topology'] == 'tail':
            mol.types[~mol.types.mask] = self.system.NonBondedTable.A2N['TAIL']
          elif mol.properties['topology'] == 'loop':
            mol.types[~mol.types.mask] = self.system.NonBondedTable.A2N['LOOP']
          elif mol.properties['topology'] == 'tie':
            mol.types[~mol.types.mask] = self.system.NonBondedTable.A2N['TIE']
  def generate_unanchored_trials(self):
    new = {}
    new['x'] = np.random.random(self.num_trials)*self.system.box.lx + self.system.box.xlo
    new['y'] = np.random.random(self.num_trials)*self.system.box.ly + self.system.box.ylo
    new['z'] = np.random.random(self.num_trials)*self.system.box.ly + self.system.box.zlo

    #since we're generating in the box, the images are all zero
    images = np.zeros(self.num_trials) 
    new['imx'] = images
    new['imy'] = images
    new['imz'] = images
    new['J']   = 1.0
    return new
  def generate_singly_anchored_trials(self,local_index,anchors,trial_xyz,trial_imxyz):
    anchor_index = anchors[local_index][0]
    if anchor_index>=self.system.nbeads:

      anchor_index -= self.system.nbeads
      if anchor_index>=local_index:
        raise ValueError('Can\'t anchor to ungrown trial bead!')

      trial_x,trial_y,trial_z = trial_xyz
      trial_imx,trial_imy,trial_imz = trial_imxyz
      anchor_x   = trial_x[0,anchor_index] #first index is arbitrary
      anchor_y   = trial_y[0,anchor_index]
      anchor_z   = trial_z[0,anchor_index]
      anchor_imx = trial_imx[0,anchor_index] #first index is arbitrary
      anchor_imy = trial_imy[0,anchor_index]
      anchor_imz = trial_imz[0,anchor_index]
    else:
      anchor_x   = self.system.x[anchor_index]
      anchor_y   = self.system.y[anchor_index]
      anchor_z   = self.system.z[anchor_index]
      anchor_imx = self.system.imx[anchor_index]
      anchor_imy = self.system.imy[anchor_index]
      anchor_imz = self.system.imz[anchor_index]

    pertubations = np.random.random((self.num_trials,3)) - 0.5
    pertubations = (pertubations.T/np.linalg.norm(pertubations,axis=1)).T

    new_x = anchor_x + pertubations[:,0]
    new_y = anchor_y + pertubations[:,1]
    new_z = anchor_z + pertubations[:,2]
    (new_x, new_y, new_z),(new_imx,new_imy,new_imz) = self.system.box.wrap_positions(new_x,new_y,new_z)

    new = {}
    new['x']   = new_x
    new['y']   = new_y
    new['z']   = new_z
    new['imx'] = new_imx + anchor_imx
    new['imy'] = new_imy + anchor_imy
    new['imz'] = new_imz + anchor_imz
    new['J']   = 1.0
    return new
  def generate_doubly_anchored_trials(self,local_index,anchors,trial_xyz,trial_imxyz):
    anchor_index1 = anchors[local_index][0]
    if anchor_index1>=self.system.nbeads:
      anchor_index1 -= self.system.nbeads
      if anchor_index1>=local_index:
        raise ValueError('Can\'t anchor to ungrown trial bead!')

      trial_x,trial_y,trial_z = trial_xyz
      trial_imx,trial_imy,trial_imz = trial_imxyz
      anchor_x1   = trial_x[0,anchor_index1] #first index is arbitrary
      anchor_y1   = trial_y[0,anchor_index1]
      anchor_z1   = trial_z[0,anchor_index1]
      anchor_imx1 = trial_imx[0,anchor_index1] #first index is arbitrary
      anchor_imy1 = trial_imy[0,anchor_index1]
      anchor_imz1 = trial_imz[0,anchor_index1]
    else:
      anchor_x1   = self.system.x[anchor_index1]
      anchor_y1   = self.system.y[anchor_index1]
      anchor_z1   = self.system.z[anchor_index1]
      anchor_imx1 = self.system.imx[anchor_index1]
      anchor_imy1 = self.system.imy[anchor_index1]
      anchor_imz1 = self.system.imz[anchor_index1]

    anchor_index2 = anchors[local_index][1]
    if anchor_index2>=self.system.nbeads:
      anchor_index2 -= self.system.nbeads
      if anchor_index2>=local_index:
        raise ValueError('Can\'t anchor to ungrown trial bead!')

      trial_x,trial_y,trial_z = trial_xyz
      trial_imx,trial_imy,trial_imz = trial_imxyz
      anchor_x2   = trial_x[0,anchor_index2] #first index is arbitrary
      anchor_y2   = trial_y[0,anchor_index2]
      anchor_z2   = trial_z[0,anchor_index2]
      anchor_imx2 = trial_imx[0,anchor_index2] #first index is arbitrary
      anchor_imy2 = trial_imy[0,anchor_index2]
      anchor_imz2 = trial_imz[0,anchor_index2]
    else:
      anchor_x2   = self.system.x[anchor_index2]
      anchor_y2   = self.system.y[anchor_index2]
      anchor_z2   = self.system.z[anchor_index2]
      anchor_imx2 = self.system.imx[anchor_index2]
      anchor_imy2 = self.system.imy[anchor_index2]
      anchor_imz2 = self.system.imz[anchor_index2]
    
    backboneVec = np.array([[anchor_x2-anchor_x1, anchor_y2-anchor_y1, anchor_z2-anchor_z1]]).T
    (new_x, new_y, new_z),_ = self.system.box.wrap_positions(backboneVec[0],backboneVec[1],backboneVec[2])
    backboneVec = [new_x[0],new_y[0],new_z[0]]
    backboneDist = np.linalg.norm(backboneVec)
    if backboneDist>=2.0:
      return None

    axis,angle = self.rotateQ.get_axis_angle(vec_from=[0.,0.,1.],vec_to=backboneVec)
    self.rotateQ.set_rotation(axis,angle)

    trial_angles = np.random.random(self.num_trials)*2*np.pi
    bond_dist = 1.0001#because floating point math is dumb
    x = np.cos(trial_angles)*(bond_dist - 0.25*backboneDist*backboneDist)**(0.5)
    y = np.sin(trial_angles)*(bond_dist - 0.25*backboneDist*backboneDist)**(0.5)
    z = [backboneDist/2.0]*self.num_trials
    baseVectors = np.array([x,y,z]).T

    new_x = np.zeros(self.num_trials)
    new_y = np.zeros(self.num_trials)
    new_z = np.zeros(self.num_trials)
    for i,vec in enumerate(baseVectors):
      x,y,z = self.rotateQ.rotate(vec)
      new_x[i] = x + anchor_x1
      new_y[i] = y + anchor_y1
      new_z[i] = z + anchor_z1
    (new_x, new_y, new_z),(new_imx,new_imy,new_imz) = self.system.box.wrap_positions(new_x,new_y,new_z)

    new = {}
    new['x']   = new_x
    new['y']   = new_y
    new['z']   = new_z
    new['imx'] = new_imx + anchor_imx1 #use anchor_imx1 because anchor_x1 is used in the loop above
    new['imy'] = new_imy + anchor_imy1
    new['imz'] = new_imz + anchor_imz1
    new['J'] = backboneDist
    return new
  def calc_rosenbluth(self,UBase,retrace=False):
    abort = False

    # trial_indices= np.arange(self.num_trials,dtype=np.int)

    if retrace:
      beacons = self.retrace_beacons
      anchors = self.retrace_anchors
    else:
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
    beacon_weights = []
    for local_index,sys_index in enumerate(self.indices):


      #################
      ## GROW_TRIALS ##
      #################
      if anchors[local_index] is None:
        new = self.generate_unanchored_trials()
      elif len(anchors[local_index])==1:
        trial_xyz = (trial_x,trial_y,trial_z)
        trial_imxyz = (trial_imx,trial_imy,trial_imz)
        new = self.generate_singly_anchored_trials(local_index,anchors,trial_xyz,trial_imxyz)
      elif len(anchors[local_index])==2:
        trial_xyz = (trial_x,trial_y,trial_z)
        trial_imxyz = (trial_imx,trial_imy,trial_imz)
        new = self.generate_doubly_anchored_trials(local_index,anchors,trial_xyz,trial_imxyz)
        if new is None:
          abort = True
          return abort,None,None,None
      else:
        raise ValueError('RosenbluthChain cannot currently handle more than 2 anchors on a growth step')

      trial_x[:,local_index]     = new['x']
      trial_y[:,local_index]     = new['y']
      trial_z[:,local_index]     = new['z']
      trial_imx[:,local_index]   = new['imx']
      trial_imy[:,local_index]   = new['imy']
      trial_imz[:,local_index]   = new['imz']
      trial_types[:,local_index] = [self.system.types[sys_index] for i in range(self.num_trials)]

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

      self.system.set_trial_move(
          x=trial_x[:,:(local_index+1)],
          y=trial_y[:,:(local_index+1)],
          z=trial_z[:,:(local_index+1)],
          imx=trial_imx[:,:(local_index+1)],
          imy=trial_imy[:,:(local_index+1)],
          imz=trial_imz[:,:(local_index+1)],
          types=trial_types[:,:(local_index+1)],
          bonds=[],
          )

      #calculation will ignore indices and individually calculate the PE for each trial bead
      U,UBond = self.engine.TPE.compute( trial_move=True, partial_indices=self.indices, ntrials=self.num_trials)
      trial_potential_energies = np.add(U,UBase)
      rosen_weights.append(np.exp(np.negative(trial_potential_energies)))

      #############
      ## BEACONS ##
      #############
      if (beacons[local_index] is not None) and (beacons[local_index][1]>1):
        end_index,nbonds = beacons[local_index]

        if end_index>=self.system.nbeads:
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

        ## Calculate guiding beacon for all trial positions
        dx = trial_x[:,local_index] - end_x
        dy = trial_y[:,local_index] - end_y
        dz = trial_z[:,local_index] - end_z
        dr = np.array([dx,dy,dz])
        all_dist = self.system.box.wrap_distances(dr[0],dr[1],dr[2])

        beaconx = self.beacon[nbonds-1]['x']
        beacony = self.beacon[nbonds-1]['y']
        beacon_values = np.interp(all_dist,beaconx,beacony,left=0,right=0)

        rosen_weights[-1] *= beacon_values
      else:
        beacon_values = np.ones(self.num_trials)


      #############
      ## SHOW IT ##
      #############
      # if not retrace and self.viz is not None:
      if self.viz is not None:
        print '===================={:02d}/{:02d}===================='.format(local_index,self.length-1)
        print 'U\n',U
        print 'UBOND\n',UBond
        print 'UBASE\n',UBase
        print 'ROSEN\n',rosen_weights[-1]
        print 'BEACONS\n',beacon_values
        #print 'CHOSEN,R,B,J',rosen_weights[-1][chosen_index],beacon_weights[-1],new['J']
        print '============================================='

        x = [self.system.x[i] for i in self.indices]
        y = [self.system.y[i] for i in self.indices]
        z = [self.system.z[i] for i in self.indices]
        orig = np.array([x,y,z]).T

        anchor_pos = []
        for i in anchors[local_index]:
          if i>=self.system.nbeads:
            i-=self.system.nbeads
            anchor_pos.append([trial_x[0,i],trial_y[0,i],trial_z[0,i]])
          else:
            anchor_pos.append([self.system.x[i],self.system.y[i],self.system.z[i]])
        anchor_pos = np.array(anchor_pos)

        if beacons[local_index] is not None:
          beacon = np.array([end_x,end_y,end_z])
        else:
          beacon = None

        x = trial_x[:,local_index]
        y = trial_y[:,local_index]
        z = trial_z[:,local_index]
        trials = np.array([x,y,z]).T

        if local_index>0:
          x = trial_x[0,:local_index]
          y = trial_y[0,:local_index]
          z = trial_z[0,:local_index]
          grown = np.array([x,y,z]).T

        else:
          grown = None

        # x = trial_x[chosen_index,local_index]
        # y = trial_y[chosen_index,local_index]
        # z = trial_z[chosen_index,local_index]
        # chosen = np.array([x,y,z])
        chosen = None

        # index_list = []
        # for mol in self.system.molecule_types['ChainSegment']:
        #   if mol.properties['topology'] == 'tail':
        #     for i in mol.properties['chain_ends']:
        #       if i not in mol.properties['connected_to']:
        #         index_list.append(i)
        # x = self.system.x[index_list]
        # y = self.system.y[index_list]
        # z = self.system.z[index_list]
        # aa = np.array([x,y,z]).T
        # anchor_pos = np.append(anchor_pos,aa,axis=0)

        self.draw_trial(anchor_pos,beacon,grown,trials,orig,chosen)



      ############
      ## CHOOSE ##
      ############
      # If all of the rosen_weights are extremely small or zero, it means that
      # the configuration is "stuck" and that all trial monomers are high energy.
      # We abort the growth early in order to not waste time continuing the growth
      # of a broken configuration.
      # if (not retrace) and np.sum(rosen_weights[-1])<1e-16:
      if (not retrace) and np.all(rosen_weights[-1]==0.):
        abort = True
        return abort,None,None,None

      if retrace:
        chosen_index = 0
      else:
        trial_probabilities = rosen_weights[-1]/np.sum(rosen_weights[-1])
        chosen_index = choice(self.num_trials,p=trial_probabilities)

      #need the chosen bias weight for the final acceptance
      beacon_weights.append(beacon_values[chosen_index])


      # Equalize all trial_values at the current step to the chosen_value
      trial_x[:,local_index]     = trial_x[chosen_index,local_index]
      trial_y[:,local_index]     = trial_y[chosen_index,local_index]
      trial_z[:,local_index]     = trial_z[chosen_index,local_index]
      trial_imx[:,local_index]   = trial_imx[chosen_index,local_index]
      trial_imy[:,local_index]   = trial_imy[chosen_index,local_index]
      trial_imz[:,local_index]   = trial_imz[chosen_index,local_index]


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

    return abort,rosen_weights,beacon_weights,new['J']
