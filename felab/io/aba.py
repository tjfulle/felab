import felab.elemlib.CPE3 as CPE3
import felab.elemlib.CPE4 as CPE4
import felab.elemlib.CPS3 as CPS3
import felab.elemlib.CPS4 as CPS4
import felab.elemlib.CAX4 as CAX4
import felab.elemlib.CPE4R as CPE4R
import felab.elemlib.CPE8 as CPE8
from felab.io.abaparse import AbaqusModel


def ElementType(name):
    name = name.upper()
    if name == "CPE4":
        return CPE4.CPE4
    elif name == "CAX4":
        return CAX4.CAX4
    elif name == "CPE4R":
        return CPE4R.CPE4R
    elif name in ("CPS4", "CPS4R"):
        return CPS4.CPS4
    elif name == "CPE3":
        return CPE3.CPE3
    elif name == "CPS3":
        return CPS3.CPS3
    if name[:4] == "CPE8":
        return CPE8.CPE8
    raise ValueError("UNKNOWN ELEMENT TYPE {0}".format(name))


def ReadInput(filename):
    model = AbaqusModel(filename)
    nodtab = model.node_table()
    eletab = model.element_table()
    element_blocks = model.element_blocks(fun=ElementType)
    surfaces = model.surfaces.todict()
    return nodtab, eletab, model.nodesets, model.elsets, surfaces, element_blocks
