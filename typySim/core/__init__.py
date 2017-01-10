from cy import CellList
from cy import Box
from cy import BondList
from Selection import *
from PairTable import *
from ValueTable import *
from System import * #must be last

# import numpy as np

# def masked_assign(masked_array,new_values):
#   if ((masked_array.mask==False).sum()) != (new_values.shape[0]):
#     raise ValueError('new_values does not match masked_array shape!')
#   value_list = list(new_values)
#   for data,mask in zip(masked_array.data,masked_array.mask):
#     if not all(mask):
#       data[:] = value_list.pop()
