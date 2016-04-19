from .abaparse import AbaqusModel

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

def ReadInput(filename):
    model = AbaqusModel(filename)
    nodtab = model.node_table()
    eletab = model.element_table()
    eleblx = model.element_blocks(fun=ElementType)
    surfaces = model.surfaces.todict()
    return nodtab, eletab, model.nodesets, model.elsets, surfaces, eleblx
