import vtk
import ipdb as pdb; st = pdb.set_trace
import numpy as np

colormap ={}
colormap[0] = (0.94,0.90,0.55)
colormap[0] = (1.0,1.0,1.0)
colormap[1] = (1.0,0.0,0.0)
colormap[2] = (1.0,0.0,1.0)
colormap[3] = (0.0,0.0,1.0)
colormap[4] = (0.0,1.0,1.0)
colormap[5] = (0.0,1.0,1.0)
colormap[6] = (0.0,0.0,0.0)
colormap = {key:(i*255,j*255,k*255) for key,(i,j,k) in colormap.items()}

class SystemViewer(object):
  def __init__(self,system):
    # create a rendering window and renderer
    self.ren = vtk.vtkRenderer()
    self.ren.SetBackground(1.0,1.0,1.0)
    self.renWin = vtk.vtkRenderWindow()
    self.renWin.AddRenderer(self.ren)
     
    # create a renderwindowinteractor
    self.iren = vtk.vtkRenderWindowInteractor()
    self.iren.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
    self.iren.SetRenderWindow(self.renWin)
    self.actor_map = {}

    self.box = None
    self.system = system
  def draw_box(self):
    if self.box is not None:
      self.ren.RemoveActor(box)

    edges = []
    edges.append(self.system.box.xlo)
    edges.append(self.system.box.xhi)
    edges.append(self.system.box.ylo)
    edges.append(self.system.box.yhi)
    edges.append(self.system.box.zlo)
    edges.append(self.system.box.zhi)

    source = vtk.vtkOutlineSource()
    source.SetBounds(*edges)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor((0.,0.,0.))
    self.ren.AddActor(actor)
    self.box = actor
  def draw_all(self):
    self.draw_bonds()
    for mol in self.system.molecules:
      self.add_molecule(mol)
  def draw_bonds(self,check_bonds=True):
    # create source
    points = vtk.vtkPoints()
    scalars= vtk.vtkUnsignedCharArray()
    scalars.SetNumberOfComponents(3)
    types = self.system.types
    pos = self.system.positions
    for t,p in zip(types,pos):
      points.InsertNextPoint(p)
      scalars.InsertNextTuple3(*colormap[t])
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.GetPointData().SetScalars(scalars)
    polyData.Allocate()

    bondArray = np.array(self.system.bonds.bonds)
    bondList = [[i,j] for i,jlist in enumerate(bondArray) for j in jlist if j!=-1]
    if bondList:
      bondList =  np.sort(bondList,axis=1)
      temp  = bondList.view(np.dtype((np.void, bondList.dtype.itemsize * bondList.shape[1])))
      _, idx = np.unique(temp, return_index=True)
      uniqueBondList = bondList[idx]
      for bond in uniqueBondList:
        b0 = bond[0]
        b1 = bond[1]
        if check_bonds==True:
          p0 = pos[b0]
          p1 = pos[b1]
          dp = np.abs(p1 - p0)
          if any(dp>self.system.box.half_L):
            continue
        b = vtk.vtkIdList()
        b.InsertNextId(bond[0])
        b.InsertNextId(bond[1])
        polyData.InsertNextCell(vtk.VTK_LINE,b)
  
      line_mapper = vtk.vtkPolyDataMapper()
      line_mapper.SetInputData(polyData)
  
      line_actor = vtk.vtkActor()
      line_actor.SetMapper(line_mapper)
      line_actor.GetProperty().SetLineWidth(3)
  
      self.ren.AddActor(line_actor)
    else:
      line_actor = None
  def draw_system(self,check_bonds=False):
    # create source
    points = vtk.vtkPoints()
    scalars= vtk.vtkUnsignedCharArray()
    scalars.SetNumberOfComponents(3)
    types = self.system.types
    pos = self.system.positions
    for t,p in zip(types,pos):
      points.InsertNextPoint(p)
      scalars.InsertNextTuple3(*colormap[t])
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.GetPointData().SetScalars(scalars)
    polyData.Allocate()


    bondList = [[i,j] for i,bonds in enumerate(self.system.bonds) for j in bonds if j!=-1]
    if bondList:
      bondList =  np.sort(bondList,axis=1)
      temp  = bondList.view(np.dtype((np.void, bondList.dtype.itemsize * bondList.shape[1])))
      _, idx = np.unique(temp, return_index=True)
      uniqueBondList = bondList[idx]
      for bond in uniqueBondList:
        b0 = bond[0]
        b1 = bond[1]
        if check_bonds==True:
          p0 = pos[b0]
          p1 = pos[b1]
          dp = np.abs(p1 - p0)
          if any(dp>self.system.box.half_L):
            continue
        b = vtk.vtkIdList()
        b.InsertNextId(bond[0])
        b.InsertNextId(bond[1])
        polyData.InsertNextCell(vtk.VTK_LINE,b)
  
      line_mapper = vtk.vtkPolyDataMapper()
      line_mapper.SetInputData(polyData)
  
      line_actor = vtk.vtkActor()
      line_actor.SetMapper(line_mapper)
      line_actor.GetProperty().SetLineWidth(3)
  
      self.ren.AddActor(line_actor)
    else:
      line_actor = None

    source = vtk.vtkSphereSource()
    beads = vtk.vtkGlyph3D()
    beads.SetColorModeToColorByScalar()
    beads.SetSourceConnection(source.GetOutputPort());
    beads.SetInputData(polyData);
    beads.ScalingOff();
    beads.Update();
     
    # mapper
    bead_mapper = vtk.vtkPolyDataMapper()
    bead_mapper.SetInputConnection(beads.GetOutputPort())
     
    # actor
    bead_actor = vtk.vtkActor()
    bead_actor.SetMapper(bead_mapper)

    # assign actor to the renderer
    self.ren.AddActor(bead_actor)

    self.actor_map['all'] = {'line':line_actor,'bead':bead_actor}
  def add_molecule(self,mol):
    # create source
    points = vtk.vtkPoints()
    scalars= vtk.vtkUnsignedCharArray()
    scalars.SetNumberOfComponents(3)
    types = mol.types.compressed()
    # pos = mol.positions.compressed().reshape(-1,3)
    x = mol.x.compressed()
    y = mol.y.compressed()
    z = mol.z.compressed()
    pos = np.array([x,y,z]).T
    for t,p in zip(types,pos):
      points.InsertNextPoint(p)
      scalars.InsertNextTuple3(*colormap[t])
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.GetPointData().SetScalars(scalars)

    source = vtk.vtkSphereSource()
    beads = vtk.vtkGlyph3D()
    beads.SetColorModeToColorByScalar()
    beads.SetSourceConnection(source.GetOutputPort());
    beads.SetInputData(polyData);
    beads.ScalingOff();
    beads.Update();
     
    # mapper
    bead_mapper = vtk.vtkPolyDataMapper()
    bead_mapper.SetInputConnection(beads.GetOutputPort())
     
    # actor
    bead_actor = vtk.vtkActor()
    bead_actor.SetMapper(bead_mapper)

    # assign actor to the renderer
    self.ren.AddActor(bead_actor)

    self.actor_map[mol] = {'bead':bead_actor}
  def show(self,picking=False):
    if picking:
      style = MouseInteractorHighLightActor()
      style.SetDefaultRenderer(self.ren)
      self.iren.SetInteractorStyle(style)
       
    self.iren.Initialize()
    self.renWin.Render()
    self.iren.Start()
