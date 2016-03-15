from numpy import *

from .utilities import *

# --------------------------------------------------------------------------- #
# -------------------------- BASE ELMENT CLASS ------------------------------ #
# --------------------------------------------------------------------------- #
class Element(object):
    nodes = None
    signature = None
    variables = None
    dimensions = None
    integration = None
    gaussp = None
    gaussw = None
    edges = []
    def jacobian(self, *args):
        raise NotImplementedError
    def shape(self, *args):
        raise NotImplementedError
    def shapegrad(self, *args):
        raise NotImplementedError
    def shapegradx(self, *args):
        raise NotImplementedError
    def response(self, *args):
        raise NotImplementedError

# --------------------------------------------------------------------------- #
# -------------------------- ELEMENT FAMILIES ------------------------------- #
# --------------------------------------------------------------------------- #
class LinknD2(Element):
    nodes = 2

class Tria3(Element):
    nodes = 3
    dimensions = 2
    edges = array([[0,1], [1,2], [2,0]])

class Tria6(Element):
    nodes = 6
    dimensions = 2
    edges = array([[0,1,3], [1,2,4], [2,0,5]])

class Quad4(Element):
    nodes = 4
    dimensions = 2
    edges = array([[0,1], [1,2], [2,3], [3,0]])

class Quad8(Element):
    nodes = 8
    dimensions = 2
    edges = array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])

def ElementFamily(dimensions, nodes, abaname=None):
    if nodes == 2:
        return LinknD2
    elif dimensions == 2 and nodes == 3:
        return Tria3
    elif dimensions == 2 and nodes == 6:
        return Tria6
    elif dimensions == 2 and nodes == 4:
        return Quad4
    elif dimensions == 2 and nodes == 8:
        return Quad8
    raise ValueError('Unknown element family')

def InterpolateToCentroid(data):
    from .isoplib2_3 import IsoPTria3
    from .isoplib2_4 import IsoPQuad4, PlaneStrainQuad4SelectiveReduced
    from .isoplib2_8 import IsoPQuad8
    numgauss, numcomp = data.shape
    if numgauss == 1:
        return data.reshape(numcomp)
    if numgauss == 3:
        return IsoPTria3.interpolate_to_centroid(data)
    if numgauss == 4:
        return IsoPQuad4.interpolate_to_centroid(data)
    if numgauss == 5:
        return PlaneStrainQuad4SelectiveReduced.interpolate_to_centroid(data)
    if numgauss == 9:
        return IsoPQuad8.interpolate_to_centroid(data)
    raise NotImplementedError
