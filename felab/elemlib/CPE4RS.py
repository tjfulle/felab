from numpy import *
from .isop_p4_base import isop_p4_base
from .isop_base import stress_displacement
from .gauss_rule_info import quad_gauss_rule_info


class CPE4RS(isop_p4_base, stress_displacement):
    ndir = 3
    nshr = 1
    num_gauss = 4
    selective_reduced = True
    srip = array([[0., 0.]])
    sriw = array([4.])

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(2, point)

    def bmatrix(self, dN, *args):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
