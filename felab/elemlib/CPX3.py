import numpy as np
from .CMDN import CMDN


# --------------------------------------------------------------------------- #
# --------------------- TRIANGLE ISOPARAMETRIC ELEMENTS --------------------- #
# --------------------------------------------------------------------------- #
class CPX3(CMDN):
    """3-node isoparametric stress-displacement element

    Notes
    -----
    Node and element face numbering


            1
            | .
       [1]  |   .  [0]
            |     .
            2------0
              [2]

    """

    nodes = 3
    dimensions = 2
    elefab = {"t": 1.0}
    cp = np.array([1.0, 1.0, 1.0]) / 3.0
    edges = np.array([[0, 1], [1, 2], [2, 0]])
    xp = np.array([[1.0, 0.0, 0], [0.0, 1.0, 0], [0.0, 0.0, 1]])
    signature = [(1, 1, 0, 0, 0, 0, 0), (1, 1, 0, 0, 0, 0, 0), (1, 1, 0, 0, 0, 0, 0)]

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        a = 0.5 * (x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]))
        return a

    @property
    def volume(self):
        return self.t * self.area

    def edge_shape(self, edge, xp):
        # ORDERING OF NODES
        xb = self.xc[self.edges[edge]]
        he = np.sqrt((xb[1, 0] - xb[0, 0]) ** 2 + (xb[1, 1] - xb[0, 1]) ** 2)
        o = np.array({0: [0, 1, 2], 1: [2, 0, 1], 2: [1, 2, 0]}[edge])
        s = he * (xp + 1) / 2.0
        return np.array([(he - s) / he, s / he, 0.0])[o]

    def shapefun_der(self, coord, qcoord):
        """Shape functions of 3 node triangle

        Parameters
        ----------
        coord : ndarray
            The coordinate in the physical coordinates
        qcoord : ndarray
            The coordinate in the triangle coordinates

        Returns
        -------
        N : ndarray
            The shape function in the natural coordinates
        dN : ndarray
            The shape function derivatives in the physical coordinates
        J : float
            The Jacobian of the transformation

        """
        x = coord[:, 0]
        y = coord[:, 1]
        z1, z2, z3 = qcoord
        # Triangle coordinates *are* the shape functions
        N = np.array([z1, z2, z3])
        dNdz = np.eye(3)
        J = np.array(
            [
                [1, 1, 1],
                [np.dot(x, dNdz[:, 0]), np.dot(x, dNdz[:, 1]), np.dot(x, dNdz[:, 2])],
                [np.dot(y, dNdz[:, 0]), np.dot(y, dNdz[:, 1]), np.dot(y, dNdz[:, 2])],
            ]
        )
        Jdet = np.linalg.det(J)
        D = np.array([[0, 0], [1, 0], [0, 1]])
        P = np.dot(np.linalg.inv(J), D)
        dNdx_T = np.dot(dNdz, P)
        return N, dNdx_T.T, Jdet
