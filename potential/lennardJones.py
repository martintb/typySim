import numpy as np


def calc(r,eps,sigma,rcut=2.5):
  if r<rcut:
    U = 4.0 * eps * ((sigma/r)**(12.0) - (sigma/r)**(6.0))
  else:
    U = 0.0
  return U
