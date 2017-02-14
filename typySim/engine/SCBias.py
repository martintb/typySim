from numpy import sign
class SCBias(object):
  def __init__(self,tail=1.0,loop=1.0,tie=1.0):
    self.tail_bias = tail
    self.loop_bias = loop
    self.tie_bias  = tie
  def factor(self,tails,loops,ties):
    factor  = 1.0
    factor *= (self.tail_bias * abs(tails))**(sign(tails))
    factor *= (self.loop_bias * abs(loops))**(sign(loops))
    factor *= (self.tie_bias  * abs(ties ))**(sign(ties ))
    return factor

