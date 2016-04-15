from numpy import *
from .isop2_8 import CSDIsoParametricQuad8 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENT --------------------- #
# --------------------------------------------------------------------------- #
c = -sqrt(3./5.)
class PlaneStrainQuad8(BaseElement):
    ndir = 3
    nshr = 1
    integration = 9
    gaussp = array([[c,  c], [0,  c], [-c,  c],
                    [c,  0], [0,  0], [-c,  0],
                    [c, -c], [0, -c], [-c, -c]])
    gaussw = array([0.30864197, 0.49382716, 0.30864197,
                    0.49382716, 0.79012346, 0.49382716,
                    0.30864197, 0.49382716, 0.30864198])
    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
