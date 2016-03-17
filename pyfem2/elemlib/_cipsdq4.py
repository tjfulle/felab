from numpy import *
from numpy.linalg import inv, det

from ._cipsd import CIPSDElement

# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENTS -------------------- #
# --------------------------------------------------------------------------- #
class CIPSDQ4Element(CIPSDElement):
    """4-node isoparametric element

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
    elefab = {'t':1.}
    signature = [(1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0)]
    dimensions = 2
    integration = 4
    gaussp = array([[-1., -1.], [ 1., -1.], [-1.,  1.], [ 1.,  1.]]) / sqrt(3.)
    gaussw = ones(4)
    cp = array([0, 0], dtype=float64)
    xp = array([[-1, -1], [1, -1], [1, 1], [-1, 1],
                [0, -1], [1, 0], [0, 1], [-1, 0]], dtype=float64),
    edges = array([[0, 1], [1, 2], [2, 3], [3, 0]])

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        A2  = (x[0]*y[1] - x[1]*y[0]) + (x[1]*y[2] - x[2]*y[1])
        A2 += (x[2]*y[3] - x[3]*y[2]) + (x[3]*y[0] - x[0]*y[3])
        return A2 / 2.

    @property
    def volume(self):
        return self.t * self.area

    def isop_map(self, xi):
        raise NotImplementedError

    def shape(self, xi, edge=None):
        if edge is not None:
            # Evaluate shape function on specific edge
            xi = array([[xi,-1.],[1.,xi],[xi,1.],[-1.,xi]][edge])
        N = array([(1. - xi[0]) * (1. - xi[1]),
                   (1. + xi[0]) * (1. - xi[1]),
                   (1. + xi[0]) * (1. + xi[1]),
                   (1. - xi[0]) * (1. + xi[1])]) / 4.
        return N

    def shapegrad(self, xi):
        dN = array([[-1. + xi[1],  1. - xi[1], 1. + xi[1], -1. - xi[1]],
                    [-1. + xi[0], -1. - xi[0], 1. + xi[0],  1. - xi[0]]]) / 4.
        return dN

    def surface_force(self, edge, qe):
        edgenod = self.edges[edge]
        xb = self.xc[edgenod]
        gw = ones(2)
        gp = array([-1./sqrt(3.), 1./sqrt(3.)])
        he = sqrt((xb[1,1]-xb[0,1])**2+(xb[1,0]-xb[0,0])**2)
        Fe = zeros(8)
        for (p, xi) in enumerate(gp):
            # Form Gauss point on specific edge
            Ne = self.shape(xi, edge=edge)
            Pe = self.pmatrix(Ne)
            Fe += he / 2. * gw[p] * dot(Pe.T, qe)
        return Fe
