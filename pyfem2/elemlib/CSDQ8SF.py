from numpy import zeros
from ._csd_q8_f import CSDQ8FElement
class PlaneStressQuad8(CSDQ8FElement):
    ndir = 2
    nshr = 1
    def bmatrix(self, dN):
        """Assemble and return the B matrix"""
        B = zeros((3, 16))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
