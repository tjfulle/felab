from numpy import *
from .CPX4 import CPX4
from .gauss_rule_info import quad_gauss_rule_info


class CPE4(CPX4):
    """4 node plane-strain stress-displacement element"""

    ndir = 3
    nshr = 1
    num_gauss = 4

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(2, point)

    def bmatrix(self, dN, *args):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
