from numpy import *
from numpy.linalg import inv, det

from ._csd_r import CSDRElement

TOOR3 = 1. / sqrt(3.)
c = -sqrt(3./5.)

class CSDQ8RElement(CSDRElement):
    """8-node isoparametric element, reduced integration

    Notes
    -----
    Node and element face numbering

              [2]
           3---6---2
           |       |
       [3] 7       5 [1]
           |       |
           0---4---1
              [0]

    """
    nodes = 8
    elefab = {'t':1.}
    signature = [(1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0),
                 (1,1,0,0,0,0,0)]
    dimensions = 2
    integration = 13
    integration1 = 9
    gaussp = array([[c,  c], [0,  c], [-c,  c],
                    [c,  0], [0,  0], [-c,  0],
                    [c, -c], [0, -c], [-c, -c],
                    [-TOOR3, -TOOR3], [ TOOR3, -TOOR3],
                    [-TOOR3,  TOOR3], [ TOOR3,  TOOR3]])
    gaussw = array([0.30864197, 0.49382716, 0.30864197,
                    0.49382716, 0.79012346, 0.49382716,
                    0.30864197, 0.49382716, 0.30864198,
                    1., 1., 1., 1.])
    cp = array([0, 0], dtype=float64)
    xp = array([[-1, -1], [1, -1], [1, 1], [-1, 1],
                [0, -1], [1, 0], [0, 1], [-1, 0]], dtype=float64),
    edges = array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])

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
        x, y = xi[:2]
        N = zeros(8)
        # CORNER NODES
        N[0] = -0.25 * (1. - x) * (1. - y) * (1. + x + y)
        N[1] =  0.25 * (1. + x) * (1. - y) * (x - y - 1.)
        N[2] =  0.25 * (1. + x) * (1. + y) * (x + y - 1.)
        N[3] =  0.25 * (1. - x) * (1. + y) * (y - x - 1.)
        # MIDSIDE NODES
        N[4] =  0.5 * (1. - x * x) * (1. - y)
        N[5] =  0.5 * (1. + x) * (1. - y * y)
        N[6] =  0.5 * (1. - x * x) * (1. + y)
        N[7] =  0.5 * (1. - x) * (1. - y * y)
        return N

    def shapegrad(self, xi):
        x, y = xi[:2]
        dN = zeros((2, 8))
        dN[0,0] =  0.25 * (1. - y) * (2. * x + y)
        dN[0,1] =  0.25 * (1. - y) * (2. * x - y)
        dN[0,2] =  0.25 * (1. + y) * (2. * x + y)
        dN[0,3] =  0.25 * (1. + y) * (2. * x - y)
        dN[0,4] = -x * (1. - y)
        dN[0,5] =  0.5 * (1. - y * y)
        dN[0,6] = -x * (1. + y)
        dN[0,7] = -0.5 * (1. - y * y)

        dN[1,0] =  0.25 * (1. - x) * (x + 2. * y)
        dN[1,1] =  0.25 * (1. + x) * (2. * y - x)
        dN[1,2] =  0.25 * (1. + x) * (2. * y + x)
        dN[1,3] =  0.25 * (1. - x) * (2. * y - x)
        dN[1,4] = -0.5 * (1. - x * x)
        dN[1,5] = -(1. + x) * y
        dN[1,6] =  0.5 * (1. - x * x)
        dN[1,7] = -(1. - x) * y
        return dN

    def surface_force(self, edge, qe):
        edgenod = self.edges[edge]
        xb = self.xc[edgenod]
        gp = array([c, 0, -c])
        gw = array([0.5555555556, 0.8888888889, 0.5555555556])
        Fe = zeros(16)
        for (p, xi) in enumerate(gp):
            # Evaluate shape function on edge
            dxdxi = dot([[-.5 + xi, .5 + xi, -2. * xi]], xb)
            Jac = sqrt(dxdxi[0, 0] ** 2 + dxdxi[0, 1] ** 2)
            Ne = self.shape(xi, edge=edge)
            Pe = self.pmatrix(Ne)
            Fe += Jac * gw[p] * dot(Pe.T, qe)
        return Fe
