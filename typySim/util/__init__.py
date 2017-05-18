import os
import copy

def write_xyz(fname,x,y,z,types):
  nbeads = len(x)
  with open(fname,'w') as f:
    f.write('{}\n'.format(nbeads))
    f.write('typySim generated XYZ file\n')
    for i in range(nbeads):
      f.write('{:>3s} {:10.5f} {:10.5f} {:10.5f}\n'.format(str(types[i]),x[i],y[i],z[i]))
    f.write('\n')
