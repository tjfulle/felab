from numpy import *
from .CPX3 import CPX3
from .gauss_rule_info import tri_gauss_rule_info


class CPS3(CPX3):
    """3 node plane-stress stress-displacement element"""

    ndir = 2
    nshr = 1
    num_gauss = 3

    @staticmethod
    def gauss_rule_info(point):
        return tri_gauss_rule_info(3, point)

    def bmatrix(self, dN, *args):
        B = zeros((3, 6))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
