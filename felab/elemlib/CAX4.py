from numpy import *
from .isop_p4_base import isop_p4_base
from .isop_base import stress_displacement
from .gauss_rule_info import quad_gauss_rule_info


class CAX4(isop_p4_base, stress_displacement):
    ndir = 3
    nshr = 1
    num_gauss = 4
    elefab = {'formulation': 1}

    @property
    def formulation(self):
        return self.axisymmetric

    @formulation.setter
    def formulation(self, arg):
        assert arg in (0, 1, 2)
        self.axisymmetric = arg

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(2, point)

    def bmatrix(self, dN, *args):
        N = args[0]
        rp = dot(N, self.xc[:,0])
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        B[2, 0::2] = N / rp
        return B
