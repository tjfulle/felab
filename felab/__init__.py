import numpy as np

from felab.mesh import *
from felab.constants import *
from felab.fe_model import FEModel
from felab.material import Material

np.set_printoptions(4, suppress=True)


def printprecision(a):
    np.set_printoptions(a, suppress=True)
