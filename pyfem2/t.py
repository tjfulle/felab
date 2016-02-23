from numpy import *
from isoplib import IsoPElement

__all__ = ['PlaneStressTria3', 'PlaneStrainTria3']

# --------------------------------------------------------------------------- #
# --------------------- TRIANGLE ISOPARAMETRIC ELEMENTS --------------------- #
# --------------------------------------------------------------------------- #
class IsoPTria3(IsoPElement):
    """3-node isoparametric element

    Notes
    -----
    Node and element face numbering


            2
            | \
       [2]  |  \ [1]
            |   \
            0----1
              [0]

    """
    signature = 1100000
    numdim, numnod, ndof = 2, 3, 2
    edges = array([[0, 1], [1, 2], [2, 0]])
    gaussp = array([[.6, .2], [.2, .6], [.2, .2]])
    gaussw = ones(3) / 6.
    cp = array([1, 1], dtype=float64) / 3.
    xp = array([[0, 0], [1, 0], [0, 1]], dtype=float64)

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
        a = .5 * (x[0] * (y[1] - y[2]) +
                  x[1] * (y[2] - y[0]) +
                  x[2] * (y[0] - y[1]))
        return a

    @property
    def volume(self):
        return self.t * self.area

    def shape(self, xi):
        Ne = array([xi[0], xi[1], 1. - xi[0] - xi[1]])
        return Ne

    def shapegrad(self, xi):
        return array([[1., 0., -1.], [0., 1., -1.]])

    def bshape(self, xi):
        return array([1. - xi, 1. + xi]) / 2.

    def bshapegrad(self, xi):
        return array([-1., 1.]) / 2.

# --------------------------------------------------------------------------- #
# ------------------------ USER ELEMENT TYPES ------------------------------- #
# --------------------------------------------------------------------------- #
class PlaneStressTria3(IsoPTria3):
    ndir, nshr = 2, 1
    def bmatrix(self, dN):
        B = zeros((3, 6))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B

class PlaneStrainTria3(IsoPTria3):
    ndir, nshr = 3, 1
    def bmatrix(self, dN):
        B = zeros((4, 6))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
