import numpy as np


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

    def __init__(self, label, elenod, elecoord, elemat):
        self.label = label
        self.inodes = np.asarray(elenod)
        self.xc = np.asarray(elecoord)
        self.material = elemat

    def eval(self, *args):
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
