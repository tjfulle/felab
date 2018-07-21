from numpy import *
from .isop_p3_base import isop_p3_base
from .isop_base import stress_displacement
from .gauss_rule_info import tri_gauss_rule_info


# --------------------------------------------------------------------------- #
# --------------------- TRIANGLE ISOPARAMETRIC ELEMENTS --------------------- #
# --------------------------------------------------------------------------- #
class CPE3(isop_p3_base, stress_displacement):
    ndir = 3
    nshr = 1
    num_gauss = 3

    @staticmethod
    def gauss_rule_info(point):
        return tri_gauss_rule_info(3, point)

    def bmatrix(self, dN, *args):
        B = zeros((4, 6))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
