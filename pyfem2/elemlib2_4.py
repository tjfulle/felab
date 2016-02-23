from numpy import *
from isoplib import IsoPElement, IsoPReduced

__all__ = ['PlaneStressQuad4', 'PlaneStrainQuad4', 'PlaneStrainQuad4Reduced',
           'PlaneStrainQuad4BBar']

# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENTS -------------------- #
# --------------------------------------------------------------------------- #
class IsoPQuad4(IsoPElement):
    """4-node isoparametric element

    Notes
    -----
    Node and element face numbering

               [2]
            3-------2
            |       |
       [3]  |       | [1]
            |       |
            0-------1
               [0]

    """
    signature = 1100000
    numdim, numnod, ndof = 2, 4, 2
    gaussp = array([[-1., -1.], [ 1., -1.], [-1.,  1.], [ 1.,  1.]]) / sqrt(3.)
    gaussw = ones(4)
    cp = array([0, 0], dtype=float64)
    xp = array([[-1, -1], [1, -1], [1, 1], [-1, 1],
                [0, -1], [1, 0], [0, 1], [-1, 0]], dtype=float64),
    bgaussw = ones(2)
    bgaussp = array([-1., 1.]) / sqrt(3.)
    edges = array([[0, 1], [1, 2], [2, 3], [3, 0]])
    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.nodes = elenod
        self.xc = elecoord
        self.material = elemat
        self.t = elefab.get('t')
        if self.t is None:
            raise ValueError('Incorrect number of element fabrication properties')

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

    def shape(self, xi):
        N = array([(1. - xi[0]) * (1. - xi[1]),
                   (1. + xi[0]) * (1. - xi[1]),
                   (1. + xi[0]) * (1. + xi[1]),
                   (1. - xi[0]) * (1. + xi[1])]) / 4.
        return N

    def shapegrad(self, xi):
        dN = array([[-1. + xi[1],  1. - xi[1], 1. + xi[1], -1. - xi[1]],
                    [-1. + xi[0], -1. - xi[0], 1. + xi[0],  1. - xi[0]]]) / 4.
        return dN

    def bshape(self, xi):
        return array([1. - xi, 1. + xi]) / 2.

    def bshapegrad(self, xi):
        return array([-1., 1.]) / 2.

# --------------------------------------------------------------------------- #
# ------------------------ USER ELEMENT TYPES ------------------------------- #
# --------------------------------------------------------------------------- #
class PlaneStressQuad4(IsoPQuad4):
    ndir, nshr = 2, 1
    def bmatrix(self, dN):
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B

class PlaneStrainQuad4(IsoPQuad4):
    ndir, nshr = 3, 1
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B

class PlaneStrainQuad4BBar(IsoPQuad4):
    ndir, nshr = 3, 1
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        # mean dilatational formulation
        dNb = self.shapegradxbar(self.xc)
        for a in range(self.numnod):
            i = 2 * a
            j = i + 1
            bb1 = (dNb[0, a] - dN[0, a]) / 2.
            bb2 = (dNb[1, a] - dN[1, a]) / 2.
            B[0, i:i+2] += [bb1, bb2]
            B[1, i:i+2] += [bb1, bb2]
        return B

class PlaneStrainQuad4Reduced(IsoPQuad4, IsoPReduced):
    ndir, nshr = 3, 1
    gaussp = array([[0., 0.]])
    gaussw = array([4.])
    hglassp = array([[0., 0.]])
    hglassv = array([[1., -1., 1., -1.]]) # hourglass vector
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
