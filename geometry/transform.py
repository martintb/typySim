import numpy as np


def normalize(v):
  v /= np.linalg.norm(v)

def rotation_matrix(v1,v2):
  R = np.empty((3,3))
  R[:,0] = normalize(v1)
  R[:,2] = normalize(np.cross(v1,v2))
  R[:,1] = normalize(np.cross(R[:,2],v1))
  return np.array(R)
  
if __name__=="__main__":
  # import vtk
  # import ipdb; ist = ipdb.set_trace
  #  
  # # create source
  # points = vtk.vtkPoints()
  # for i in range(10):
  #   points.InsertNextPoint([0.,0.,i*1.0])
  # polyData = vtk.vtkPolyData()
  # polyData.SetPoints(points)

  # source = vtk.vtkSphereSource()

  # glyph3D = vtk.vtkGlyph3D()
  # glyph3D.SetSourceConnection(source.GetOutputPort());
  # glyph3D.SetInputData(polyData);
  # glyph3D.Update();
  #  
  # # mapper
  # mapper = vtk.vtkPolyDataMapper()
  # mapper.SetInputConnection(glyph3D.GetOutputPort())
  #  
  # # actor
  # actor = vtk.vtkActor()
  # actor.SetMapper(mapper)
  # actor.GetProperty().SetColor(1.0,0.0,0.0)
  #  
  # # create a rendering window and renderer
  # ren = vtk.vtkRenderer()
  # ren.SetBackground(1.0,1.0,1.0)
  # renWin = vtk.vtkRenderWindow()
  # renWin.AddRenderer(ren)
  #  
  # # create a renderwindowinteractor
  # iren = vtk.vtkRenderWindowInteractor()
  # iren.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
  # iren.SetRenderWindow(renWin)

  # # assign actor to the renderer
  # ren.AddActor(actor)
  #  
  # # enable user interface interactor
  # iren.Initialize()
  # renWin.Render()
  # iren.Start()
  from ..viz import Scene
  scene = Scene()
  pos1 = [[0.,0.,i] for i in range(10)]
  pos2 = [[0.,i,0.] for i in range(10)]
  scene.add_beads(pos1)
  scene.add_beads(pos2)
  scene.show()

