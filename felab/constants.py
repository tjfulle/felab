import os
from numpy import sqrt, array, float64, zeros
# MAXIMUM NUMBER OF POSSIBLE DOF PER NODE
MDOF = 7

# FREEDOM LABELS
X, Y, Z = 0, 1, 2  # DISPLACEMENTS
TX, TY, TZ = 3, 4, 5  # ROTATIONS
T = 6  # TEMPERATURE

Rr, Zr = 0, 1  # DISPLACEMENTS

# BOUNDARY LABELS
ALL = '_All'
PIN = '_Pin'
ILO, IHI = '_Ilo', '_Ihi'
JLO, JHI = '_Jlo', '_Jhi'
KLO, KHI = '_Klo', '_Khi'
BOUNDARY = '_Boundary'
NEUMANN = 'Neumann'
DIRICHLET = 'Dirichlet'

# RADIAL
RLO, RHI = '_Ilo', '_Ihi'
ZLO, ZHI = '_Jlo', '_Jhi'

# EDGES
S1, S2, S3, S4, S5, S6, S7, S8, S9, S10 = range(10)

TOL1 = 1e-12

# FIELD TYPES
SCALAR = 'Scalar'
VECTOR = 'Vector'
SYMTENSOR = 'Symmetric Tensor'
TENSOR = 'Tensor'

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

# TENSORS
ONED = '1D'
TWOD = '3D'
THREED = '3D'
PSTRAIN = 'Plane Strain'
PSTRESS = 'Plane Stress'

# DISTRIBUTED LOADS
SLOAD = 1
DLOAD = 2
SFLUX = 3
SFILM = 4
HSRC = 5

SMALL_DISPLACEMENT = 0
LARGE_DISPLACEMENT = 1

# PROCEDURE TYPES
STATIC_ITERATIVE = 1
STATIC_DIRECT = 2
DYNAMIC = 17
HEAT_TRANSFER_STEADY_STATE = 31

STIFF_AND_RHS = 1
STIFF_ONLY = 2
MASS_ONLY = 4
RHS_ONLY = 5
MASS_AND_RHS = 6

GENERAL = 0
LINEAR_PERTURBATION = 1
