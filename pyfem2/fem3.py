from numpy import *
from numpy.linalg import solve, LinAlgError

from .constants import *
from .utilities import *
from .fem1 import FiniteElementModel
from .elemlib2_3T import DiffussiveHeatTransfer2D3

class HeatTransfer2DModel(FiniteElementModel):
    numdim = 2
    def init(self):
        # Request allocation of field variables
        self.request_output_variable('T', SCALAR, NODE)
        self.request_output_variable('R', SCALAR, NODE)

    def Solve(self):
        self.validate(DiffussiveHeatTransfer2D3)
        K = self.assemble_global_stiffness(self.sfilm)
        F, Q = self.assemble_global_force(self.src, self.sflux, self.sfilm)
        Kbc, Fbc = self.apply_bc(K, F+Q)
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')
        Ft = dot(K, self.dofs)
        R = Ft - F - Q
        self.snapshot(T=self.dofs, R=R)
        self.T = self.dofs
