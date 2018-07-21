from numpy import *
from .isop_p8_base import isop_p8_base
from .isop_base import stress_displacement
from .gauss_rule_info import quad_gauss_rule_info


class CPS8(isop_p8_base, stress_displacement):
    ndir = 2
    nshr = 1
    num_gauss = 9

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(rule=3, point=point)

    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = zeros((3, 16))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
