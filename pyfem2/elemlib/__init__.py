from .heat_transfer import *
from .stress_displacement import *
from .link import *
from .beam import *

# --------------------------------------------------------------------------- #
# -------------------------- ELEMENT FAMILIES ------------------------------- #
# --------------------------------------------------------------------------- #
from .element import Element as _Element
from numpy import array as _array
class LinknD2(_Element):
    nodes = 2
    edges = _array([])

class Tria3(_Element):
    nodes = 3
    dimensions = 2
    edges = _array([[0,1], [1,2], [2,0]])

class Tria6(_Element):
    nodes = 6
    dimensions = 2
    edges = _array([[0,1,3], [1,2,4], [2,0,5]])

class Quad4(_Element):
    nodes = 4
    dimensions = 2
    edges = _array([[0,1], [1,2], [2,3], [3,0]])

class Quad8(_Element):
    nodes = 8
    dimensions = 2
    edges = _array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])

def ElementFamily(dimensions, nodes, ndir=None, nshr=None, abaname=None):
    if nodes == 2:
        return LinknD2
    elif dimensions == 2 and nodes == 3:
        return PlaneStressTria3  # Tria3
    elif dimensions == 2 and nodes == 6:
        return Tria6
    elif dimensions == 2 and nodes == 4:
        return PlaneStressQuad4  # Quad4
    elif dimensions == 2 and nodes == 8:
        return PlaneStressQuad8  # Quad8
    raise ValueError('Unknown element family')
