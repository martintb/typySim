import numpy as np


def calc(r,eps,sigma,rcut=2.0**(1.0/6.0)):
  if r<rcut:
    U = 4.0 * eps * ((sigma/r)**(12.0) - (sigma/r)**(6.0)) + eps
  else:
    U = 0.0
  return U
