from numpy import *
from numpy.linalg import det, inv
from .continuum_stress_disp_quad8_full import CSDQ8FElement
class PlaneStrainQuad8BBar(CSDQ8FElement):
    ndir = 3
    nshr = 1
    def bmatrix(self, dN):
        """Assemble and return the B matrix"""
        B = zeros((4, 16))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]

        # MEAN DILATATIONAL FORMULATION
        xc = self.xc

        dNb = zeros((2,8))
        jac = zeros(self.integration)
        for p in range(self.integration):
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
            B[0, i:i+2] += [bb1, bb2]
            B[1, i:i+2] += [bb1, bb2]

        return B
