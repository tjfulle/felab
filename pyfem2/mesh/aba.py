from .abaparse import AbaqusModel
from ..material import Material
from ..utilities import UserInputError

def ElementType(name):
    import pyfem2.elemlib as elemlib
    name = name.upper()
    if name == 'CPE4':
        return elemlib.PlaneStrainQuad4
    elif name == 'CAX4':
        return elemlib.AxiSymmetricQuad4
    elif name == 'CPE4R':
        return elemlib.PlaneStrainQuad4Reduced
    elif name in ('CPS4', 'CPS4R'):
        return elemlib.PlaneStressQuad4
    elif name == 'CPE3':
        return elemlib.PlaneStrainTria3
    elif name == 'CPS3':
        return elemlib.PlaneStressTria3
    if name[:4] == 'CPE8':
        return elemlib.PlaneStrainQuad8
    raise ValueError('UNKNOWN ELEMENT TYPE {0}'.format(name))

def abamat_to_pyfem2(dict_):
    """Parse the abaqus material definition and convert to pyfem2"""
    name = dict_['name']
    kwds = {}
    for (key, d) in dict_.items():
        if key in ('name', 'data'):
            continue
        elif key == 'elastic':
            data = [float(x) for x in d['data'][0]]
            kwds['elastic'] = dict(E=data[0], Nu=data[1])
        elif key == 'expansion':
            kwds['expansion'] = float(d['data'][0])
        elif key == 'density':
            kwds['density'] = float(d['data'][0])
        else:
            raise UserInputError('UNKNOWN ABAQUS MATERIAL KEYWORD {0}'.format(key))
    return Material(name, **kwds)

def ReadInput(filename):
    model = AbaqusModel(filename)
    nodtab = model.node_table()
    eletab = model.element_table()
    eleblx = model.element_blocks(ele_fun=ElementType, mat_fun=abamat_to_pyfem2)
    surfaces = model.surfaces.todict()
    return nodtab, eletab, model.nodesets, model.elsets, surfaces, eleblx
