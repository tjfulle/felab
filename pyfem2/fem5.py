""".pft.py: Plane Beam Column Truss application code.

"""
from numpy import *
from numpy.linalg import solve, LinAlgError

from .constants import *
from .utilities import linsolve
from .fem1 import FiniteElementModel
from .elemlib import ElasticLink2D2, PlaneBeamColumn

# --------------------------------------------------------------------------- #
# -------------------------- APPLICATION CODE ------------------------------- #
# --------------------------------------------------------------------------- #
class PlaneBeamColumnTrussModel(FiniteElementModel):
    dimensions = 2

    def Solve(self):
        # active DOF set dynamically
        self.setup((ElasticLink2D2, BeamColumn2D))

        # Assemble the global stiffness and force
        flags = [1, 0, 1, 1]
        du = zeros(self.numdof)
        K, rhs = self.assemble(self.dofs, du)
        Kbc, Fbc = self.apply_bc(K, rhs)

        # Solve
        self.dofs[:] = linsolve(Kbc, Fbc)

        # Total force, including reaction, and reaction
        R = dot(K, self.dofs) - rhs

        Q = self.external_force_array()
        R = zeros(self.numdof)
        U, R, Urot, Rrot = self.format_displacements_and_reactions(self.dofs, R)

        frame = self.steps.last.Frame(1.)
        frame.field_outputs['U'].add_data(U)
        frame.field_outputs['R'].add_data(R)
        frame.converged = True

        self.u = U

    def format_displacements_and_reactions(self, u, R):
        # construct displacement and rotation vectors
        ut, rt = zeros((self.numnod, 2)), zeros((self.numnod, 2))
        ur, rr = zeros(self.numnod), zeros(self.numnod)
        for n in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[n,j] > 0:
                    ii = self.nodfmt[n] + ix
                    if j <= 2:
                        ut[n,j] = u[ii]
                        rt[n,j] = R[ii]
                    else:
                        # This would need modifying if 3D beams were allowed
                        ur[n] = u[ii]
                        rr[n] = R[ii]
                    ix += 1
        return ut, rt, ur, rr
