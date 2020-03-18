import numpy as np
from .CMDN import CMDN


# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENTS -------------------- #
# --------------------------------------------------------------------------- #
class CPX4(CMDN):
    """4-node isoparametric element stress-displacement element base

    Notes
    -----
    Node and element face numbering

              [2]
           3-------2
           |       |
       [3] |       | [1]
           |       |
           0-------1
              [0]

    """

    nodes = 4
    dimensions = 2
    elefab = {"t": 1.0}
    signature = [
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
    ]
    cp = np.array([0, 0], dtype=float)
    xp = np.array(
        [[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0]],
        dtype=float,
    )
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        A2 = (x[0] * y[1] - x[1] * y[0]) + (x[1] * y[2] - x[2] * y[1])
        A2 += (x[2] * y[3] - x[3] * y[2]) + (x[3] * y[0] - x[0] * y[3])
        return A2 / 2.0

    @property
    def volume(self):
        return self.t * self.area

    def shape(self, qcoord, edge=None):
        if edge is not None:
            # EVALUATE SHAPE FUNCTION ON SPECIFIC EDGE
            xi = qcoord
            qcoord = np.array([[xi, -1.0], [1.0, xi], [xi, 1.0], [-1.0, xi]][edge])
        xi, eta = qcoord
        N = (
            np.array(
                [
                    (1.0 - xi) * (1.0 - eta),
                    (1.0 + xi) * (1.0 - eta),
                    (1.0 + xi) * (1.0 + eta),
                    (1.0 - xi) * (1.0 + eta),
                ]
            )
            / 4.0
        )
        return N

    def shapegrad(self, qcoord):
        xi, eta = qcoord
        dNdxi = (
            np.array(
                [
                    [-1.0 + eta, 1.0 - eta, 1.0 + eta, -1.0 - eta],
                    [-1.0 + xi, -1.0 - xi, 1.0 + xi, 1.0 - xi],
                ]
            )
            / 4.0
        )
        return dNdxi

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
