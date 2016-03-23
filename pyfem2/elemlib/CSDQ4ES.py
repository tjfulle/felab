from numpy import *
from .continuum_stress_disp_quad4_sel_red import CSDQ4SElement

class PlaneStrainQuad4SelectiveReduced(CSDQ4SElement):
    ndir = 3
    nshr = 1
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
