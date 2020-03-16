import numpy as np
from .CPX3 import CPX3
from .gauss_rule_info import tri_gauss_rule_info


class CPE3(CPX3):
    """3 node plane-strain stress-displacement element"""

    ndir = 3
    nshr = 1
    num_gauss = 3

    @staticmethod
    def gauss_rule_info(point):
        return tri_gauss_rule_info(3, point)

    def bmatrix(self, dN, *args):
        B = np.zeros((4, 6))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
