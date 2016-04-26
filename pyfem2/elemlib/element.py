from numpy import *
from ..utilities import *

# --------------------------------------------------------------------------- #
# ------------------------- BASE ELEMENT CLASS ------------------------------ #
# --------------------------------------------------------------------------- #
class Element(object):
    ndir = None
    nshr = None
    nodes = None
    elefab = None
    signature = None
    dimensions = None
    integration = None
    edges = []

    def __init__(self, label, elenod, elecoord, elemat, **elefab):

        self.label = label
        self.inodes = asarray(elenod)
        self.xc = asarray(elecoord)
        self.material = elemat

        if self.elefab is None:
            if elefab:
                raise UserInputError('Element takes no element '
                                     'fabrication properties')
        else:
            unknown = [key for key in elefab if key not in self.elefab]
            if unknown:
                raise UserInputError('Unrecognized element fabrication '
                                     'properties: {0}'.format(','.join(unknown)))
            for (name, default) in self.elefab.items():
                p = elefab.get(name, default)
                if p is None:
                    raise UserInputError('Missing required fabrication '
                                         'property {0}'.format(name))
                setattr(self, name, p)

    def response(self, *args):
        raise NotImplementedError

    @classmethod
    def interpolate_to_centroid(cls, *args):
        raise NotImplementedError

    @classmethod
    def variables(cls):
        return None

    def varidx(self, key):
        return [x[0] for x in self.variables()].index(key)

    def numalloc(self, disp=1):
        m = 0
        for v in self.variables():
            name, vtype = v[:2]
            if vtype == SCALAR:
                m += 1
            elif vtype == VECTOR:
                m += model.dimensions
            elif vtype == SYMTENSOR:
                m += self.ndir + self.nshr
            elif vtype == TENSOR:
                m += self.ndir + 2*self.nshr
            else:
                raise Exception('UNKNOWN VARIABLE TYPE')
        if disp and self.integration:
            m *= self.integration
        return m
