from numpy import *
from .isoplib import CSDIsoParametricElement as BaseElement

# --------------------------------------------------------------------------- #
# --------------------- TRIANGLE ISOPARAMETRIC ELEMENTS --------------------- #
# --------------------------------------------------------------------------- #
class CSDIsoParametricTria3(BaseElement):
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
    elefab = {'t':1.}
    cp = array([1., 1.]) / 3.
    edges = array([[0, 1], [1, 2], [2, 0]])
    xp = array([[1., 0.], [0., 1.], [0., 0.]])
    signature = [(1,1,0,0,0,0,0), (1,1,0,0,0,0,0), (1,1,0,0,0,0,0)]

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        a = .5 * (x[0] * (y[1] - y[2]) +
                  x[1] * (y[2] - y[0]) +
                  x[2] * (y[0] - y[1]))
        return a

    @property
    def volume(self):
        return self.t * self.area

    def shape(self, xi, edge=None):
        if edge is not None:
            # EVALUATE SHAPE FUNCTION ON EDGE
            if   edge == 0: xi = [1.-xi/sqrt(2.), xi/sqrt(2.)]
            elif edge == 1: xi = [0., xi]
            elif edge == 2: xi = [xi, 0.]
        Ne = array([xi[0], xi[1], 1. - xi[0] - xi[1]])
        return Ne

    def shapegrad(self, xi):
        return array([[1., 0., -1.], [0., 1., -1.]])
