from numpy import *
from numpy.linalg import inv, det
from .isop_p4_base import isop_p4_base
from .isop_base import stress_displacement
from .gauss_rule_info import quad_gauss_rule_info


class CPS4I(isop_p4_base, stress_displacement):
    ndir = 2
    nshr = 1
    incompatible_modes = True
    num_gauss = 4

    @staticmethod
    def gauss_rule_info(self, point=None):
        return quad_gauss_rule_info(2, point)

    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B

    def gmatrix(self, xi):
        """Assemble and return the G matrix"""
        # ALGORITHM IN
        # THE FINITE ELEMENT METHOD: ITS BASIS AND FUNDAMENTALS
        # BY OLEK C ZIENKIEWICZ, ROBERT L TAYLOR, J.Z. ZHU

        # JACOBIAN AT ELEMENT CENTROID
        # COMPUTE THE SHAPE FUNCTION AT THE CENTROID
        dN0dxi = self.shapegrad(self.cp)

        # COMPUTE THE DEFORMATION GRADIENT AT CENTROID
        # AND THE JACOBIAN
        dx0dxi = dot(dN0dxi, self.xc)
        dxidx0 = inv(dx0dxi)
        J0 = det(dx0dxi)

        # COMPUTE THE JACOBIAN OF THE ELEMENT
        _, _, J = self.shapefun_der(xi)

        # COMPUTE DNDXI ASSOCIATED WITH THE INCOMPATIBLE MODES AND THEN FROM IT
        # AND THE JACOBIANS COMPUTED ABOVE COMPUTE DNDX
        # N = [1 - xi**2, 1 - eta**2]
        dNdxi = array([[-2*xi[0], 0], [0, -2*xi[1]]])
        dNdx = J0 / J * dot(dxidx0, dNdxi)

        # FORM THE INCOMPATIBLE G MATRIX
        G = array([[dNdx[0,0], 0,         dNdx[0,1], 0],
                   [0,         dNdx[1,0], 0,         dNdx[1,1]],
                   [dNdx[0,1], dNdx[0,0], dNdx[1,1], dNdx[1,0]]])

        return G
