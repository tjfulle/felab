from .abaparse import AbaqusModel
from .. import elemlib


def ElementType(name):
    name = name.upper()
    if name == "CPE4":
        return elemlib.CPE4
    elif name == "CAX4":
        return elemlib.CAX4
    elif name == "CPE4R":
        return elemlib.CPE4R
    elif name in ("CPS4", "CPS4R"):
        return elemlib.CPS4
    elif name == "CPE3":
        return elemlib.CPE3
    elif name == "CPS3":
        return elemlib.CPS3
    if name[:4] == "CPE8":
        return elemlib.CPE8
    raise ValueError("UNKNOWN ELEMENT TYPE {0}".format(name))


def ReadInput(filename):
    model = AbaqusModel(filename)
    nodtab = model.node_table()
    eletab = model.element_table()
    element_blocks = model.element_blocks(fun=ElementType)
    surfaces = model.surfaces.todict()
    return nodtab, eletab, model.nodesets, model.elsets, surfaces, element_blocks
