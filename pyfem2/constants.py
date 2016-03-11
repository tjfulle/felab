import os
from numpy import sqrt, array, float64, zeros
# Maximum number of possible DOF per node
MDOF = 7

# FREEDOM LABELS
X, Y, Z = 0, 1, 2  # DISPLACEMENTS
TX, TY, TZ = 3, 4, 5  # ROTATIONS
T = 6  # TEMPERATURE

# BOUNDARY LABELS
ALL = 'ALL'
ILO, IHI = 'ILO', 'IHI'
JLO, JHI = 'JLO', 'JHI'
KLO, KHI = 'KLO', 'KHI'
BOUNDARY = 'BOUNDARY'
NEUMANN = 0
DIRICHLET = 1

# EDGES
S1, S2, S3, S4, S5, S6, S7, S8, S9, S10 = range(10)

TOL1 = 1e-12

# FIELD TYPES
SCALAR = 'Scalar'
VECTOR = 'Vector'
SYMTENSOR = 'Tensor'

# SOLVERS
NEWTON = 'Newton'
RIKS = 'Riks'

# FIELD POSITIONS
NODE = 'Node'
ELEMENT = 'Element'
ELEMENT_CENTROID = 'Element Centroid'
INTEGRATION_POINT = 'Integration Point'

# MISC. CONSTANT ARRAYS AND SCALARS
Z6 = zeros(6)
I6 = array([1., 1., 1., 0., 0., 0.])
I9 = array([1., 0., 0., 0., 1., 0., 0., 0., 1.])
VOIGHT = array([1, 1, 1, 2, 2, 2], dtype=float64)

DEFAULT_TEMP = 298.

ROOT2 = sqrt(2.0)
ROOT3 = sqrt(3.0)
TOOR2 = 1.0 / ROOT2
TOOR3 = 1.0 / ROOT3
ROOT23 = ROOT2 / ROOT3

_d = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_D = os.path.join(_d, 'data')
if not os.path.isdir(DATA_D):
    import logging
    logging.warn('pyfem2 data directory not located')
