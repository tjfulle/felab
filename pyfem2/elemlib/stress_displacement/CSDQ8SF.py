from numpy import *
from .isop2_8 import CSDIsoParametricQuad8 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# --------------------------------------------------------------------------- #
c = -sqrt(3./5.)
class PlaneStressQuad8(BaseElement):
    ndir = 2
    nshr = 1
    integration = 9
    gaussp = array([[c,  c], [0,  c], [-c,  c],
                    [c,  0], [0,  0], [-c,  0],
                    [c, -c], [0, -c], [-c, -c]])
    gaussw = array([0.30864197, 0.49382716, 0.30864197,
                    0.49382716, 0.79012346, 0.49382716,
                    0.30864197, 0.49382716, 0.30864198])
    def bmatrix(self, dN):
        """Assemble and return the B matrix"""
        B = zeros((3, 16))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
