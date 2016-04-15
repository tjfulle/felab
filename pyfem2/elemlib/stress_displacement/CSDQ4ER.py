from numpy import *
from .isop2_4 import CSDIsoParametricQuad4 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# -------------------------- REDUCED INTEGRATION ---------------------------- #
# --------------------------------------------------------------------------- #
class PlaneStrainQuad4Reduced(BaseElement):
    ndir = 3
    nshr = 1
    integration = 1
    gaussw = array([4.])
    gaussp = array([[0., 0.]])
    hourglass_control = True
    # HOURGLASS CONTROL
    hglassp = array([[0., 0.]])
    hglassv = array([[1., -1., 1., -1.]])
    def bmatrix(self, dN, *args):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
