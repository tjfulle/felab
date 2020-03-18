import os
import re
import felab.util.tty as tty


# Patterns to ignore in the commands directory when looking for commands.
ignore_files = (
    r"^\.|^__init__.py$|^#|flycheck_|^_.*\.py$|"
    r"^element_base.py$|^gauss_rule_info.py$|^isop_base.py$"
)


#: global, cached list of all elements -- access through all_elements()
_all_elements = None


def all_elements():
    """Get a sorted list of all felab elements.

    It does not actually import the python files -- just gets the names.
    """
    global _all_elements
    if _all_elements is None:
        _all_elements = []
        element_paths = [os.path.dirname(os.path.realpath(__file__))]
        for path in element_paths:
            for file in os.listdir(path):
                if file.endswith(".py") and not re.search(ignore_files, file):
                    element = re.sub(r".py$", "", file)
                    _all_elements.append(element)
        _all_elements.sort()
    return _all_elements


def get_module(elem_name):
    """Imports the module for a particular element and returns it.

    Parameters
    ----------
    elem_name : str
        name of the element for which to get a module
    """
    # Import the element
    module_name = "{0}.{1}".format(__name__, elem_name)
    module = __import__(module_name, fromlist=[elem_name], level=0)
    if not hasattr(module, elem_name):  # pragma: no cover
        tty.die(
            "Element module {0} ({1}) must define class {2!r}.".format(
                module.__name__, module.__file__, elem_name
            )
        )
    return module


def Element(*, type):
    for element in all_elements():
        if type.upper() == element.upper():
            module = get_module(element)
            return getattr(module, element)
    tty.die(f"Unknown element type {type}")


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
        return Tria3
    elif dimensions == 2 and nodes == 6:
        return Tria6
    elif dimensions == 2 and nodes == 4:
        return Quad4
    elif dimensions == 2 and nodes == 8:
        return Quad8
    tty.die("Unknown element family")
