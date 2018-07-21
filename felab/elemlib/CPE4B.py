from numpy import *
from numpy.linalg import det, inv
from .isop_p4_base import isop_p4_base
from .isop_base import stress_displacement
from .gauss_rule_info import quad_gauss_rule_info


# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# ---------------------- MEAN DILATATIONAL FORUMULA ------------------------- #
# --------------------------------------------------------------------------- #
class CPE4B(isop_p4_base, stress_displacement):
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

        # MEAN DILATATIONAL FORMULATION
        xc = self.xc

        dNb = zeros((2,4))
        jac = zeros(self.num_gauss)
        for p in range(self.num_gauss):
            # COMPUTE THE INTEGRALS OVER THE VOLUME
            xi = self.gaussp[p]
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            jac[p] = det(dxdxi)
            dxidx = inv(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            dNb += dNdx * self.gaussw[p] * jac[p] / self.dimensions

        # SCALING
        ev = dot(jac, self.gaussw)
        dNb /= ev

        for a in range(self.nodes):
            i = 2 * a
            j = i + 1
            bb1 = (dNb[0, a] - dN[0, a]) / 2.
            bb2 = (dNb[1, a] - dN[1, a]) / 2.
            B[0,i:i+2] += [bb1, bb2]
            B[1,i:i+2] += [bb1, bb2]

        return B
