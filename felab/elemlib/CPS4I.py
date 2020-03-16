import numpy as np
from numpy.linalg import inv, det
from .CPX4 import CPX4
from .gauss_rule_info import quad_gauss_rule_info


class CPS4I(CPX4):
    """4 node plane-stress stress-displacement element with incompatible modes"""

    ndir = 2
    nshr = 1
    incompatible_modes = True
    num_gauss = 4

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(2, point)

    def bmatrix(self, dN, *args):
        """Assemble and return the B matrix"""
        B = np.zeros((3, 8))
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
        dx0dxi = np.dot(dN0dxi, self.xc)
        dxidx0 = inv(dx0dxi)
        J0 = det(dx0dxi)

        # COMPUTE THE JACOBIAN OF THE ELEMENT
        _, _, J = self.shapefun_der(self.xc, xi)

        # COMPUTE DNDXI ASSOCIATED WITH THE INCOMPATIBLE MODES AND THEN FROM IT
        # AND THE JACOBIANS COMPUTED ABOVE COMPUTE DNDX
        # N = [1 - xi**2, 1 - eta**2]
        dNdxi = np.array([[-2 * xi[0], 0], [0, -2 * xi[1]]])
        dNdx = J0 / J * np.dot(dxidx0, dNdxi)

        # FORM THE INCOMPATIBLE G MATRIX
        G = np.array(
            [
                [dNdx[0, 0], 0, dNdx[0, 1], 0],
                [0, dNdx[1, 0], 0, dNdx[1, 1]],
                [dNdx[0, 1], dNdx[0, 0], dNdx[1, 1], dNdx[1, 0]],
            ]
        )

        return G
