from numpy import *

# --------------------------------------------------------------------------- #
# ------------------------- BASE ELEMENT CLASS ------------------------------ #
# --------------------------------------------------------------------------- #
class Element(object):
    elefab = {}
    nodes = None
    variables = None
    signature = None
    dimensions = None
    integration = None
    edges = []

    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.inodes = asarray(elenod)
        self.xc = asarray(elecoord)
        self.material = elemat
        unknown = [key for key in elefab if key not in self.elefab]
        if unknown:
            raise UserInputError('Unrecognized element fabrication '
                                 'properties: {0}'.format(','.join(unknown)))
        for (name, default) in self.elefab.items():
            p = elefab.get(name, default)
            setattr(self, name, p)

    def response(self, *args):
        raise NotImplementedError

    @classmethod
    def interpolate_to_centroid(cls, *args):
        raise NotImplementedError
