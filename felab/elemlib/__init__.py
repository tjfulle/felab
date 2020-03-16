from .DC2D3 import DC2D3

DiffussiveHeatTransfer2D3 = DC2D3

from .CPE3 import CPE3
from .CPS3 import CPS3

from .CPE4B import CPE4B
from .CPE4 import CPE4
from .CPE4R import CPE4R
from .CPE4RS import CPE4RS

from .CAX4 import CAX4

from .CPS4 import CPS4
from .CPS4I import CPS4I

from .CPE8B import CPE8B
from .CPE8 import CPE8
from .CPE8R import CPE8R

from .CPS8 import CPS8

from .L1D2 import L1D2
from .L2D2 import L2D2
from .L3D2 import L3D2

from .B2D2 import B2D2

# --------------------------------------------------------------------------- #
# -------------------------- ELEMENT FAMILIES ------------------------------- #
# --------------------------------------------------------------------------- #
from .element_base import element_base as _Element
from numpy import array as _array


class LinknD2(_Element):
    nodes = 2
    edges = _array([])


class Tria3(_Element):
    nodes = 3
    dimensions = 2
    edges = _array([[0, 1], [1, 2], [2, 0]])


class Tria6(_Element):
    nodes = 6
    dimensions = 2
    edges = _array([[0, 1, 3], [1, 2, 4], [2, 0, 5]])


class Quad4(_Element):
    nodes = 4
    dimensions = 2
    edges = _array([[0, 1], [1, 2], [2, 3], [3, 0]])


class Quad8(_Element):
    nodes = 8
    dimensions = 2
    edges = _array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])


def element_family(dimensions, nodes, ndir=None, nshr=None, abaname=None):
    if nodes == 2:
        return LinknD2
    elif dimensions == 2 and nodes == 3:
        return CPS3  # Tria3
    elif dimensions == 2 and nodes == 6:
        return Tria6
    elif dimensions == 2 and nodes == 4:
        return CPS4  # Quad4
    elif dimensions == 2 and nodes == 8:
        return CPS8  # Quad8
    raise ValueError("Unknown element family")
