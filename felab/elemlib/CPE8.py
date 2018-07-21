from numpy import *
from .isop_p8_base import isop_p8_base
from .isop_base import stress_displacement
from .gauss_rule_info import quad_gauss_rule_info


class CPE8(isop_p8_base, stress_displacement):
    ndir = 3
    nshr = 1
    num_gauss = 9

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(3, point)

    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
