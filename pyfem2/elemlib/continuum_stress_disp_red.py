import logging
from numpy import *
from numpy.linalg import det, inv

from .continuum_stress_disp import CSDElement
from ..utilities import *

# --------------------------------------------------------------------------- #
# -------------- REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------------- #
# --------------------------------------------------------------------------- #
class CSDRElement(CSDElement):
    """Continue isoparametric stress displacement element, reduced integration"""
    gaussp = None
    gaussw = None
    variables = ('E', 'DE', 'S')
    integration = None
    hourglass_control = False

    def response(self, u, du, time, dtime, kstep, kframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type):
        """Assemble the element stiffness"""

        response = super(CSDFElement, self)._response(*args)

        if not self.hourglass_control:
            return response

        if cflag not in (STIFF_AND_FORCE, STIFF_ONLY):
            return response

        if cflag == STIFF_AND_FORCE:
            Ke = response[0]
            rhs = response[1]

        elif cflag == STIFF_ONLY:
            Ke = response

        # PERFORM HOURGLASS CORRECTION
        xc = self.xc  # + u.reshape(self.xc.shape)
        Khg = zeros(Ke.shape)
        for p in range(len(self.hglassp)):

            # SHAPE FUNCTION DERIVATIVE AT HOURGLASS GAUSS POINTS
            xi = array(self.hglassp[p])
            dNdxi = self.shapegrad(xi)

            # JACOBIAN TO NATURAL COORDINATES
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            J = det(dxdxi)

            # HOURGLASS BASE VECTORS
            g = array(self.hglassv[p])
            for i in range(len(xi)):
                xi[i] = dot(g, xc[:,i])

            # CORRECT THE BASE VECTORS TO ENSURE ORTHOGONALITY
            scale = 0.
            for a in range(self.nodes):
                for i in range(self.dimensions):
                    g[a] -= xi[i] * dNdx[i,a]
                    scale += dNdx[i,a] * dNdx[i,a]
            scale *= .01 * self.material.G

            for a in range(self.nodes):
                n1 = count_digits(self.signature[a])
                for i in range(n1):
                    for b in range(self.nodes):
                        n2 = count_digits(self.signature[b])
                        for j in range(n2):
                            K = n1 * a + i
                            L = n2 * b + j
                            Khg[K,L] += scale * g[a] * g[b] * J * 4.

        Ke += Khg

        if cflag == STIFF_ONLY:
            return Ke

        elif cflag == STIFF_AND_FORCE:
            return Ke, rhs
