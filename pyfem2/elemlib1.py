from numpy import *

from .utilities import *

# --------------------------------------------------------------------------- #
# -------------------------- BASE ELMENT CLASS ------------------------------ #
# --------------------------------------------------------------------------- #
class Element(object):
    nfab = 0
    ndof, numnod, numdim = None, None, None
    gaussp, gaussw = None, None
    signature = None
    edges = []
    def jacobian(self, *args):
        raise NotImplementedError
    def shape(self, *args):
        raise NotImplementedError
    def shapegrad(self, *args):
        raise NotImplementedError
    def shapegradx(self, *args):
        raise NotImplementedError
    def stiffness(self, *args):
        raise NotImplementedError
    def force(self, *args):
        raise NotImplementedError

# --------------------------------------------------------------------------- #
# -------------------------- ELEMENT FAMILIES ------------------------------- #
# --------------------------------------------------------------------------- #
class LinknD2(Element):
    numnod = 2

class Tria3(Element):
    numdim, numnod = 2, 3
    edges = array([[0,1], [1,2], [2,0]])

class Tria6(Element):
    numdim, numnod = 2, 6
    edges = array([[0,1,3], [1,2,4], [2,0,5]])

class Quad4(Element):
    numdim, numnod = 2, 4
    edges = array([[0,1], [1,2], [2,3], [3,0]])

class Quad8(Element):
    numdim, numnod = 2, 8
    edges = array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])

def ElementFamily(numdim, numnod, abaname=None):
    if numnod == 2:
        return LinknD2
    elif numdim == 2 and numnod == 3:
        return Tria3
    elif numdim == 2 and numnod == 6:
        return Tria6
    elif numdim == 2 and numnod == 4:
        return Quad4
    elif numdim == 2 and numnod == 8:
        return Quad8
    raise ValueError('Unknown element family')

def InterpolateToCentroid(data):
    from .isoplib2_3 import IsoPTria3
    from .isoplib2_4 import IsoPQuad4
    from .isoplib2_8 import IsoPQuad8
    numgauss, numcomp = data.shape
    if numgauss == 1:
        return data.reshape(numcomp)
    if numgauss == 3:
        return IsoPTria3.interpolate_to_centroid(data)
    if numgauss == 4:
        return IsoPQuad4.interpolate_to_centroid(data)
    if numgauss == 9:
        return IsoPQuad8.interpolate_to_centroid(data)
    raise NotImplementedError
