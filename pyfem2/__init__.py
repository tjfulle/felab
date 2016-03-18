import os
import sys
import argparse
import logging
from numpy import *
from .constants import *
from .elemlib import *
from .mesh.exodusii import File, PutNodalSolution
from .first_fe_program import UniformBar
from .finite_element_model import FiniteElementModel
from .fem2 import TrussModel
from .fem5 import PlaneBeamColumnTrussModel
from .mat import Material
from .mesh.mesh import Mesh

set_printoptions(4, suppress=True)
def printprecision(a):
    set_printoptions(a, suppress=True)

p = argparse.ArgumentParser()
p.add_argument('-v', default=1, type=int,
               help='Level of verbosity [default: %(default)s]')
p.add_argument('-p', default=4, type=int,
               help='Print precision [default: %(default)s')
args, other = p.parse_known_args(sys.argv[1:])

loglevel = {-2: logging.CRITICAL,
            -1: logging.ERROR,
             0: logging.WARNING,
             1: logging.INFO,
             2: logging.DEBUG}.get(max(min(args.v,2),-2))
logging.basicConfig(level=loglevel, format='%(levelname)s:%(message)s')

printprecision(args.p)
