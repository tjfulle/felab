from numpy import *
from .continuum_stress_disp_tria3_full import CSDT3FElement
class PlaneStressTria3(CSDT3FElement):
    ndir = 2
    nshr = 1
    def bmatrix(self, dN):
        B = zeros((3, 6))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
