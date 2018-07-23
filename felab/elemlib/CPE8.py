from numpy import *
from .CPX8 import CPX8
from .gauss_rule_info import quad_gauss_rule_info


class CPE8(CPX8):
    """8 node plane-strain element"""
    ndir = 3
    nshr = 1
    num_gauss = 9

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(rule=3, point=point)

    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
