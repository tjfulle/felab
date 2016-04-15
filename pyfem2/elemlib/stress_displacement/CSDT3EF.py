from numpy import *
from .isop2_3 import CSDIsoParametricTria3 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- TRIANGLE ISOPARAMETRIC ELEMENTS --------------------- #
# --------------------------------------------------------------------------- #
class PlaneStrainTria3(BaseElement):
    ndir = 3
    nshr = 1
    integration = 3
    gaussp = array([[.6, .2], [.2, .6], [.2, .2]])
    gaussw = ones(3) / 6.
    def bmatrix(self, dN, *args):
        B = zeros((4, 6))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
