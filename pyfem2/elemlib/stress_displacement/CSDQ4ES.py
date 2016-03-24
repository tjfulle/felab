from numpy import *
from .isop2_4 import CSDIsoParametricQuad4 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# --------------------- SELECTIVE REDUCED INTEGRATION ----------------------- #
# --------------------------------------------------------------------------- #
class PlaneStrainQuad4SelectiveReduced(BaseElement):
    ndir = 3
    nshr = 1
    integration = 4
    gaussw = ones(4)
    gaussp = array([[-1., -1.], [ 1., -1.], [-1.,  1.], [ 1.,  1.]]) / sqrt(3.)
    selective_reduced = True
    srip = array([[0., 0.]])
    sriw = array([4.])
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
