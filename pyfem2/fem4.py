from numpy import *
from numpy.linalg import solve, LinAlgError, det

from .constants import *
from .utilities import *
from .fem1 import FiniteElementModel
from .isoplib import IsoPElement

class Plane2DModel(FiniteElementModel):
    numdim = 2
    def init(self):
        self.request_output_variable('U', VECTOR, NODE)
        self.request_output_variable('R', VECTOR, NODE)
        self.request_output_variable('S', SYMTENSOR, INTEGRATION_POINT)
        self.request_output_variable('E', SYMTENSOR, INTEGRATION_POINT)
        self.request_output_variable('DE', SYMTENSOR, INTEGRATION_POINT)

    def Solve(self):
        self.validate(IsoPElement)
        K = self.assemble_global_stiffness()
        F, Q = self.assemble_global_force(self.dload, self.sload)
        Kbc, Fbc = self.apply_bc(K, F+Q)
        try:
            u = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')
        Ft = dot(K, u)
        self.dofs[:] = u.reshape(self.dofs.shape)
        R = Ft - F - Q
        R = R.reshape(self.mesh.coord.shape)
        u = u.reshape(self.mesh.coord.shape)
        self.u = u
        self.update_state()
        self.snapshot(U=u, R=R)

    def update_state(self):
        # compute the element stiffness and scatter to global array
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            E = self.steps[-1].field_outputs['E'][ieb]
            S = self.steps[-1].field_outputs['S'][ieb]
            for (e, xel) in enumerate(eb.labels):
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                u = self.dofs[el.nodes]
                de1, e1, s1 = el.update_state(u, E[e], S[e])
                self.steps[-1].field_outputs['DE'].data[ieb][iel] = de1
                self.steps[-1].field_outputs['E'].data[ieb][iel] = e1
                self.steps[-1].field_outputs['S'].data[ieb][iel] = s1
        return
