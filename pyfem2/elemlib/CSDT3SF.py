from numpy import zeros
from .csd_t3_f import CSDT3FElement
class PlaneStressTria3(CSDT3FElement):
    ndir = 2
    nshr = 1
    def bmatrix(self, dN):
        B = zeros((3, 6))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
