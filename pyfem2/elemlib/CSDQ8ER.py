from numpy import *
from ._csd_q8_r import CSDQ8RElement

TOOR3 = 1. / sqrt(3.)

class PlaneStrainQuad8Reduced(CSDQ8RElement):
    ndir = 3
    nshr = 1
    def bmatrix(self, dN):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
