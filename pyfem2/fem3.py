from numpy import *
from numpy.linalg import solve, LinAlgError

from .constants import *
from .utilities import *
from .fem1 import FiniteElementModel
from .elemlib2_3T import DiffussiveHeatTransfer2D3

class HeatTransfer2DModel(FiniteElementModel):
    dimensions = 2
    def Solve(self):
        self.setup(DiffussiveHeatTransfer2D3)
        flags = [21, 0, 1, 1]
        du = zeros(self.numdof)
        K, rhs = self.assemble(self.dofs, du, procedure=HEAT_TRANSFER)
        Kbc, Fbc = self.apply_bc(K, rhs)
        self.dofs[:] = linsolve(Kbc, Fbc)
        Q = dot(K, self.dofs) - rhs

        # Create new frame to hold updated state
        frame = self.steps.last.Frame(1.)
        frame.field_outputs['T'].add_data(self.dofs)
        frame.field_outputs['Q'].add_data(Q)
        frame.converged = True
        self.T = self.dofs
