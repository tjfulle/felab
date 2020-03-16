import numpy as np
from .CPX4 import CPX4
from .gauss_rule_info import quad_gauss_rule_info


class CPE4R(CPX4):
    """4 node plane-strain stress-displacement element with reduced integration"""

    ndir = 3
    nshr = 1
    num_gauss = 1
    elefab = {"t": 1.0, "hourglass_control": True}

    # HOURGLASS CONTROL
    hglassp = np.array([[0.0, 0.0]])
    hglassv = np.array([[1.0, -1.0, 1.0, -1.0]])

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(1, point)

    def bmatrix(self, dN, *args):
        B = np.zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
