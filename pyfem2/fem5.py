""".pft.py: Plane Beam Column Truss application code.

"""
from numpy import *
from numpy.linalg import solve, LinAlgError

from .constants import *
from .fem1 import FiniteElementModel
from .elemlibN_2 import ElasticLink2D2, BeamColumn2D

# --------------------------------------------------------------------------- #
# -------------------------- APPLICATION CODE ------------------------------- #
# --------------------------------------------------------------------------- #
class PlaneBeamColumnTrussModel(FiniteElementModel):
    numdim = 2

    def Solve(self):
        # active DOF set dynamically
        self.validate((ElasticLink2D2, BeamColumn2D))

        # Assemble the global stiffness and force
        K = self.assemble_global_stiffness()
        F, Q = self.assemble_global_force()
        Kbc, Fbc = self.apply_bc(K, F+Q)

        # Solve
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')

        # Total force, including reaction, and reaction
        Ft = dot(K, self.dofs)
        R = Ft - F - Q

        U, R, Urot, Rrot = self.format_displacements_and_reactions(self.dofs, R)

        frame = self.steps.last.Frame(1.)
        frame.field_outputs['U'].add_data(U)
        frame.field_outputs['R'].add_data(R)

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
