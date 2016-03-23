from numpy import *
from .continuum_stress_disp_tria3_full import CSDT3FElement
class PlaneStrainTria3(CSDT3FElement):
    ndir = 3
    nshr = 1
    def bmatrix(self, dN):
        B = zeros((4, 6))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
