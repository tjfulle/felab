import logging
from numpy import *
from numpy.linalg import det, inv

from .continuum_stress_disp import CSDElement
from ..utilities import *

# --------------------------------------------------------------------------- #
# ----------------FULLY INTEGRATED ISOPARAMETRIC ELEMENTS ------------------- #
# --------------------------------------------------------------------------- #
class CSDFElement(CSDElement):
    """Continuum Isoparametric Stress Displacement Element"""
    gaussp = None
    gaussw = None
    variables = ('E', 'DE', 'S')
    integration = None

    def response(self, *args):
        """Assemble the element stiffness"""
        return super(CSDFElement, self)._response(*args)
