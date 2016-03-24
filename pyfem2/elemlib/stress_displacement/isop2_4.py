from numpy import *
from .isoplib import CSDIsoParametricElement as BaseElement

# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENTS -------------------- #
# --------------------------------------------------------------------------- #
class CSDIsoParametricQuad4(BaseElement):
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
    dimensions = 2
    elefab = {'t':1.}
    signature = [(1,1,0,0,0,0,0), (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0), (1,1,0,0,0,0,0)]
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

    def shape(self, xi, edge=None):
        if edge is not None:
            # EVALUATE SHAPE FUNCTION ON SPECIFIC EDGE
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
