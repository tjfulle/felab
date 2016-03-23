from numpy import *
from numpy.linalg import inv, det

from .continuum_stress_disp_full import CSDFElement

# --------------------------------------------------------------------------- #
# --------------------- TRIANGLE ISOPARAMETRIC ELEMENTS --------------------- #
# --------------------------------------------------------------------------- #
class CSDT3FElement(CSDFElement):
    """3-node isoparametric element

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
    elefab = {'t':1.}
    signature = [(1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0)]
    dimensions = 2
    integration = 3
    gaussp = array([[.6, .2], [.2, .6], [.2, .2]])
    gaussw = ones(3) / 6.
    cp = array([1., 1.]) / 3.
    xp = array([[1., 0.], [0., 1.], [0., 0.]])
    edges = array([[0, 1], [1, 2], [2, 0]])

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

    def surface_force(self, edge, qe):
        edgenod = self.edges[edge]
        xb = self.xc[edgenod]
        gw = ones(2)
        gp = array([-1./sqrt(3.), 1./sqrt(3.)])
        he = sqrt((xb[1,1]-xb[0,1])**2+(xb[1,0]-xb[0,0])**2)
        Fe = zeros(6)
        for (p, xi) in enumerate(gp):
            # EVALUATE SHAPE FUNCTION ON EDGE. SINCE EDGES ARE NOT DEFINED ON
            # [-1, 1], THE GAUSS POINT MUST BE TRANSFORMED.
            Ne = self.shape(he/2.+xi*he/2., edge=edge)
            Pe = self.pmatrix(Ne)
            Fe += he / 2. * gw[p] * dot(Pe.T, qe)
        return Fe
