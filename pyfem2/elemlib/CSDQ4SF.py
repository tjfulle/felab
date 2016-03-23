from numpy import *

from .continuum_stress_disp_quad4_full import CSDQ4FElement

# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRESS ELEMENT ----------------------- #
# --------------------------------------------------------------------------- #
class PlaneStressQuad4(CSDQ4FElement):
    ndir = 2
    nshr = 1
    def bmatrix(self, dN):
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
