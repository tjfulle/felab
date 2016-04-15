from numpy import *
from .isop2_8 import CSDIsoParametricQuad8 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# -------------------------- REDUCED INTEGRATION ---------------------------- #
# --------------------------------------------------------------------------- #
class PlaneStrainQuad8Reduced(BaseElement):
    ndir = 3
    nshr = 1
    integration = 4
    gaussp = array([[-sqrt(3.), -sqrt(3.)], [ sqrt(3.), -sqrt(3.)],
                    [-sqrt(3.),  sqrt(3.)], [ sqrt(3.),  sqrt(3.)]])
    gaussw = array([1., 1., 1., 1.])
    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
