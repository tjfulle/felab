from numpy import *
from .isop2_4 import CSDIsoParametricQuad4 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRESS ELEMENT ----------------------- #
# --------------------------------------------------------------------------- #
class PlaneStressQuad4(BaseElement):
    ndir = 2
    nshr = 1
    integration = 4
    gaussw = ones(4)
    gaussp = array([[-1., -1.], [ 1., -1.], [-1.,  1.], [ 1.,  1.]]) / sqrt(3.)
    def bmatrix(self, dN):
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
