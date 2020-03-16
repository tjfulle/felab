import numpy as np
from numpy.linalg import det, inv
from .CPX4 import CPX4
from .gauss_rule_info import quad_gauss_rule_info


class CPE4B(CPX4):
    """4 node plane-strain stress-displacement element with bbar stabilization"""

    ndir = 3
    nshr = 1
    num_gauss = 4

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(2, point)

    def bmatrix(self, dN, *args):
        B = np.zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]

        # MEAN DILATATIONAL FORMULATION
        dNb = np.zeros((2, 4))
        w = np.zeros(self.num_gauss)
        jac = np.zeros(self.num_gauss)
        for p in range(self.num_gauss):
            # COMPUTE THE INTEGRALS OVER THE VOLUME
            xi, w[p] = self.gauss_rule_info(p)
            _, dN, jac[p] = self.shapefun_der(self.xc, xi)
            dNb += dN * w[p] * jac[p] / self.dimensions

        # SCALING
        ev = np.dot(jac, w)
        dNb /= ev

        for a in range(self.nodes):
            i = 2 * a
            bb1 = (dNb[0, a] - dN[0, a]) / 2.0
            bb2 = (dNb[1, a] - dN[1, a]) / 2.0
            B[0, i : i + 2] += [bb1, bb2]
            B[1, i : i + 2] += [bb1, bb2]

        return B
