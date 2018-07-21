from numpy import *
from ..utilities import *

# --------------------------------------------------------------------------- #
# ------------------------- BASE ELEMENT CLASS ------------------------------ #
# --------------------------------------------------------------------------- #
class element_base(object):
    ndir = None
    nshr = None
    nodes = None
    elefab = None
    signature = None
    dimensions = None
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

    @staticmethod
    def num_integration():
        return None

    @classmethod
    def interpolate_to_centroid(cls, *args):
        raise NotImplementedError

    @staticmethod
    def variables():
        return None
