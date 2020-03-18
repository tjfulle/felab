import numpy as np
from .CMDN import CMDN


# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENTS -------------------- #
# --------------------------------------------------------------------------- #
class CPX8(CMDN):
    """8 node isoparametric stress-displacement element

    Notes
    -----
    Node and element face numbering

               [2]
            3---6---2
            |       |
       [3]  7       5 [1]
            |       |
            0---4---1
               [0]

    """

    nodes = 8
    dimensions = 2
    elefab = {"t": 1.0}
    cp = np.array([0, 0], dtype=float)
    xp = np.array(
        [[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0]],
        dtype=float,
    )
    edges = np.array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])
    signature = [
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
    ]

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        A2 = (x[0] * y[1] - x[1] * y[0]) + (x[1] * y[2] - x[2] * y[1])
        A2 += (x[2] * y[3] - x[3] * y[2]) + (x[3] * y[0] - x[0] * y[3])
        return A2 / 2.0

    @property
    def volume(self):
        return self.t * self.area

    # def foo():
    #     if edge is not None:
    #         # EVALUATE SHAPE FUNCTION ON SPECIFIC EDGE
    #         xi = np.array([[xi, -1.0], [1.0, xi], [xi, 1.0], [-1.0, xi]][edge])

    def shape(self, qcoord, edge=None):
        if edge is not None:
            # EVALUATE SHAPE FUNCTION ON SPECIFIC EDGE
            x = qcoord
            qcoord = np.array([[x, -1.0], [1.0, x], [x, 1.0], [-1.0, x]][edge])
        xi, eta = qcoord[:2]
        # CORNER NODES
        N = np.zeros(8)
        N[0] = -0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 + xi + eta)
        N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0)
        N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0)
        N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) * (eta - xi - 1.0)
        # MIDSIDE NODES
        N[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta)
        N[5] = 0.5 * (1.0 + xi) * (1.0 - eta * eta)
        N[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta)
        N[7] = 0.5 * (1.0 - xi) * (1.0 - eta * eta)
        return N

    def shapegrad(self, qcoord):
        xi, eta = qcoord[:2]
        dN = np.zeros((2, 8))
        dN[0, 0] = 0.25 * (1.0 - eta) * (2.0 * xi + eta)
        dN[0, 1] = 0.25 * (1.0 - eta) * (2.0 * xi - eta)
        dN[0, 2] = 0.25 * (1.0 + eta) * (2.0 * xi + eta)
        dN[0, 3] = 0.25 * (1.0 + eta) * (2.0 * xi - eta)
        dN[0, 4] = -xi * (1.0 - eta)
        dN[0, 5] = 0.5 * (1.0 - eta * eta)
        dN[0, 6] = -xi * (1.0 + eta)
        dN[0, 7] = -0.5 * (1.0 - eta * eta)

        dN[1, 0] = 0.25 * (1.0 - xi) * (xi + 2.0 * eta)
        dN[1, 1] = 0.25 * (1.0 + xi) * (2.0 * eta - xi)
        dN[1, 2] = 0.25 * (1.0 + xi) * (2.0 * eta + xi)
        dN[1, 3] = 0.25 * (1.0 - xi) * (2.0 * eta - xi)
        dN[1, 4] = -0.5 * (1.0 - xi * xi)
        dN[1, 5] = -(1.0 + xi) * eta
        dN[1, 6] = 0.5 * (1.0 - xi * xi)
        dN[1, 7] = -(1.0 - xi) * eta
        return dN

    def shapefun_der(self, coord, qcoord):
        """Shape function and derivative of 4 node bilinear element

        Parameters
        ----------
        coord : ndarray
            The coordinate in the physical coordinates
        qcoord : ndarray
            The coordinate in the natural coordinates

        Returns
        -------
        N : ndarray
            The shape function evaluated at Gauss coordinate
        dNdx : ndarray
            The shape function derivative in the physical coordinates
        J : float
            The Jacobian of the transformation

        """

        # Shape function and derivative at Gauss point
        N = self.shape(qcoord)
        dNdxi = self.shapegrad(qcoord)

        # Jacobian to natural coordinates
        dxdxi = np.dot(dNdxi, coord)
        dxidx = np.linalg.inv(dxdxi)
        J = np.linalg.det(dxdxi)

        # Convert shape function derivatives to derivatives wrt global physical
        # coordinates
        dNdx = np.dot(dxidx, dNdxi)

        return N, dNdx, J
