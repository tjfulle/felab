import os
import sys
from numpy import *
from .constants import *
from .elemlib import *
from .exodusii import File, PutNodalSolution
from .fem0 import UniformBar
from .fem1 import FiniteElementModel
from .fem2 import TrussModel
from .fem3 import HeatTransfer2DModel
from .fem4 import Plane2DModel
from .fem5 import PlaneBeamColumnTrussModel
from .mat import Material
from .mesh import Mesh

set_printoptions(1, suppress=True)
def printprecision(a):
    set_printoptions(a)
