import vtk

class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):
  def __init__(self,parent=None):
    self.AddObserver("LeftButtonPressEvent",self.leftButtonPressEvent)
    self.LastPickedActor = None
    self.LastPickedProperty = vtk.vtkProperty()
  def leftButtonPressEvent(self,obj,event):
    clickPos = self.GetInteractor().GetEventPosition()
 
    picker = vtk.vtkPropPicker()
    picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
 
    # get the new
    self.NewPickedActor = picker.GetActor()
 
    # If something was selected
    if self.NewPickedActor:
      # If we picked something before, reset its property
      if self.LastPickedActor:
          self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)
 
 
      # Save the property of the picked actor so that we can
      # restore it next time
      self.LastPickedProperty.DeepCopy(self.NewPickedActor.GetProperty())
      # Highlight the picked actor by changing its properties
      self.NewPickedActor.GetProperty().SetColor(1.0, 0.0, 0.0)
      self.NewPickedActor.GetProperty().SetDiffuse(1.0)
      self.NewPickedActor.GetProperty().SetSpecular(0.0)
 
      # save the last picked actor
      self.LastPickedActor = self.NewPickedActor
    self.OnLeftButtonDown()
    return

class MouseInteractorHighLightActor2(vtk.vtkInteractorStyleTrackballCamera):
  def __init__(self,ren=None,parent=None):
    self.AddObserver("LeftButtonPressEvent",self.leftButtonPressEvent)
    self.LastActorCopy = None
    if ren is not None:
      self.ren=ren
  def leftButtonPressEvent(self,obj,event):
    clickPos = self.GetInteractor().GetEventPosition()
 
    picker = vtk.vtkPropPicker()
    picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
 
    # get the new
    self.NewPickedActor = picker.GetActor()
 
    # If something was selected
    if self.NewPickedActor:
      self.ActorCopy = vtk.vtkActor().Deepcopy(self.NewPickedActor)
      if self.LastActorCopy:
          self.LastActorCopy.Delete()
 
      # Highlight the picked actor by changing its properties
      self.ActorCopy.GetProperty().SetRepresentationToWireFrame()
      self.ActorCopy.GetProperty().SetColor(1.0, 0.0, 0.0)
      self.ActorCopy.GetProperty().SetDiffuse(1.0)
      self.ActorCopy.GetProperty().SetSpecular(0.0)
 
      # save the last picked actor
      self.LastActorCopy = self.ActorCopy
    self.OnLeftButtonDown()
    return
