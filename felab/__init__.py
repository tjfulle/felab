import numpy as np

from felab.elemlib import *

from felab.constants import *
from felab.first_fe_program import UniformBar
from felab.fe_model import fe_model
from felab.material import Material
from felab.mesh import *
from felab.mesh.exodusii import File, put_nodal_solution

np.set_printoptions(4, suppress=True)


def printprecision(a):
    np.set_printoptions(a, suppress=True)
