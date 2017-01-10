from math import sqrt
import numpy as np

def line(length,diameter):
  molData = {}
  molData['positions'] = [[0.,0.,i*diameter] for i in range(length)]
  # molData['bonds'] = [[i,j] for i,j in zip(range(length-1),range(1,length))]
  # molData['angles'] = [[i,j,k] for i,j,k in zip(range(length-2),range(1,length-1),range(2,length))]
  return molData
