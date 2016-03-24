from numpy import *
from .isop2_3 import CSDIsoParametricTria3 as BaseElement
class PlaneStressTria3(BaseElement):
    ndir = 2
    nshr = 1
    integration = 3
    gaussp = array([[.6, .2], [.2, .6], [.2, .2]])
    gaussw = ones(3) / 6.
    def bmatrix(self, dN):
        B = zeros((3, 6))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
