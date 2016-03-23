import logging
from numpy import *
from numpy.linalg import det, inv

from .continuum_stress_disp import CSDElement
from ..utilities import *

# --------------------------------------------------------------------------- #
# ---------------- INCOMPATIBLE MODES ISOPARAMETRIC ELEMENTS ---------------- #
# --------------------------------------------------------------------------- #
class CSDIElement(CSDElement):
    """Continuum isoparametric stress displacement element, incompatible modes"""
    incompatible_modes = True

    def response(self, *args):
        return super(CSDIElement, self)._response(*args)
