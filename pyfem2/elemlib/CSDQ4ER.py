from numpy import *
from ._csd_q4_r import CSDQ4RElement

class PlaneStrainQuad4Reduced(CSDQ4RElement):
    ndir = 3
    nshr = 1
    # HOURGLASS CONTROL
    hglassp = array([[0., 0.]])
    hglassv = array([[1., -1., 1., -1.]])
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
