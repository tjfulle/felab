from numpy import *
from .continuum_stress_disp_quad8_full import CSDQ8FElement
class PlaneStrainQuad8(CSDQ8FElement):
    ndir = 3
    nshr = 1
    def bmatrix(self, dN):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
