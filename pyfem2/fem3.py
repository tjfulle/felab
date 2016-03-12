from numpy import *
from numpy.linalg import solve, LinAlgError

from .constants import *
from .utilities import *
from .fem1 import FiniteElementModel
from .elemlib2_3T import DiffussiveHeatTransfer2D3

class HeatTransfer2DModel(FiniteElementModel):
    numdim = 2
    def Solve(self):
        self.setup(DiffussiveHeatTransfer2D3)
        K = self.assemble_global_stiffness(self.sfilm)
        F, Q = self.assemble_global_force(self.src, self.sflux, self.sfilm)
        Kbc, Fbc = self.apply_bc(K, F+Q)
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')
        Ft = dot(K, self.dofs)
        R = Ft - F - Q

        # Create new frame to hold updated state
        frame = self.steps.last.Frame(1.)
        frame.field_outputs['T'].add_data(self.dofs)
        frame.field_outputs['R'].add_data(R)
        frame.converged = True
        self.T = self.dofs
