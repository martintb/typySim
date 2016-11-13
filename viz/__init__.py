import vtk
# import numpy as np
import ipdb; ist = ipdb.set_trace

class MolecularViewer(object):
  colormap ={}
  colormap[0] = (1.0,0.0,0.0)
  colormap[1] = (0.0,0.5,0.0)
  colormap[2] = (0.0,0.0,1.0)
  colormap[3] = (1.0,1.0,0.0)
  colormap[4] = (1.0,0.0,1.0)
  colormap[5] = (0.0,1.0,1.0)
  colormap[6] = (0.0,0.0,0.0)
  colormap[7] = (1.0,1.0,1.0)
  colormap = {key:(i*255,j*255,k*255) for key,(i,j,k) in colormap.items()}
  def __init__(self,L=None):
    # create a rendering window and renderer
    self.ren = vtk.vtkRenderer()
    self.ren.SetBackground(1.0,1.0,1.0)
    self.renWin = vtk.vtkRenderWindow()
    self.renWin.AddRenderer(self.ren)
     
    # create a renderwindowinteractor
    self.iren = vtk.vtkRenderWindowInteractor()
    self.iren.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
    self.iren.SetRenderWindow(self.renWin)
    self.groups = []

    if L is not None:
      L = [-L/2.0 if i%2==0 else +L/2.0 for i in range(6)]
      source = vtk.vtkOutlineSource()
      source.SetBounds(*L)
      mapper = vtk.vtkPolyDataMapper()
      mapper.SetInputConnection(source.GetOutputPort())
      actor = vtk.vtkActor()
      actor.SetMapper(mapper)
      actor.GetProperty().SetColor((0.,0.,0.))
      self.ren.AddActor(actor)

      self.box = actor
    else:
      self.box = None
  def add_beads(self,pos,color=(0.0,0.0,1.0)):
    # create source
    points = vtk.vtkPoints()
    for p in pos:
      points.InsertNextPoint(p)
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    source = vtk.vtkSphereSource()

    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetSourceConnection(source.GetOutputPort());
    glyph3D.SetInputData(polyData);
    glyph3D.Update();
     
    # mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph3D.GetOutputPort())
     
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)

    # assign actor to the renderer
    self.ren.AddActor(actor)

    group = {}
    group['actor'] = actor
    group['mapper'] = mapper
    group['glyph3D'] = glyph3D
    group['polyData'] = polyData
    group['points'] = points
    self.groups.append(group)
  def add_molecule(self,mol):
    # create source
    points = vtk.vtkPoints()
    scalars= vtk.vtkUnsignedCharArray()
    scalars.SetNumberOfComponents(3)
    types = mol.types.compressed()
    pos = mol.positions.compressed().reshape(-1,3)
    for t,p in zip(types,pos):
      points.InsertNextPoint(p)
      scalars.InsertNextTuple3(*self.colormap[t])
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.GetPointData().SetScalars(scalars)

    source = vtk.vtkSphereSource()
    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetColorModeToColorByScalar()
    glyph3D.SetSourceConnection(source.GetOutputPort());
    glyph3D.SetInputData(polyData);
    glyph3D.ScalingOff();
    glyph3D.Update();
     
    # mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph3D.GetOutputPort())
     
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # assign actor to the renderer
    self.ren.AddActor(actor)

    group = {}
    group['actor'] = actor
    group['mapper'] = mapper
    group['glyph3D'] = glyph3D
    group['polyData'] = polyData
    group['points'] = points
    self.groups.append(group)
  def show(self,picking=False):
    if picking:
      style = MouseInteractorHighLightActor()
      style.SetDefaultRenderer(self.ren)
      self.iren.SetInteractorStyle(style)
       
    self.iren.Initialize()
    self.renWin.Render()
    self.iren.Start()
  def translate_test_1(self):
    import ipdb; ist = ipdb.set_trace
    import time
               
    transform = vtk.vtkTransform()
    transform.PostMultiply()
    transform.Translate(1.0, 0, 0);
    
    # enable user interface interactor
    for i in range(10):
      transform.Translate(1.0, 0, 0);
      self.groups[0]['actor'].SetUserTransform(transform)
      self.renWin.Render()
      time.sleep(0.1)

    self.iren.Initialize()
    self.iren.Start()
  def translate_test_2(self):
    import ipdb; ist = ipdb.set_trace
    import time

    points = self.groups[0]['points']
               
    # enable user interface interactor
    for i in range(10):
      for j in range(points.GetNumberOfPoints()):
        p1 = points.GetPoint(j)
        p = [k+1.0 for k in p1]
        print p1,p
        points.SetPoint(j,p)
      points.Modified()
      self.renWin.Render()
      time.sleep(0.1)

    self.iren.Initialize()
    self.iren.Start()
