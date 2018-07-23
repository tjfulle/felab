from numpy import *
from numpy.linalg import det, inv
from .CPX8 import CPX8
from .gauss_rule_info import quad_gauss_rule_info


class CPE8B(CPX8):
    """8 node plane-strain element with bbar stabilization"""
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

        # MEAN DILATATIONAL FORMULATION
        xc = self.xc

        dNb = zeros((2,8))
        w = zeros(self.num_gauss)
        jac = zeros(self.num_gauss)
        for p in range(self.num_gauss):
            # COMPUTE THE INTEGRALS OVER THE VOLUME
            xi, w[p] = self.gauss_rule_info(p)
            _, dN, jac[p] = self.shapefun_der(self.xc, xi)
            dNb += dN * w[p] * jac[p] / self.dimensions

        # SCALING
        ev = dot(jac, w)
        dNb /= ev

        for a in range(self.nodes):
            i = 2 * a
            j = i + 1
            bb1 = (dNb[0, a] - dN[0, a]) / 2.
            bb2 = (dNb[1, a] - dN[1, a]) / 2.
            B[0, i:i+2] += [bb1, bb2]
            B[1, i:i+2] += [bb1, bb2]

        return B
