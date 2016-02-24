from numpy import array, dot, zeros
from numpy.linalg import solve, LinAlgError

from .fem1 import FiniteElementModel
from .elemlibN_2 import ElasticLinknD2
from .constants import *

__all__ = ['TrussModel']

class TrussModel(FiniteElementModel):
    numdim = None

    def init(self):
        # Request allocation of field variables
        self.request_output_variable('U', VECTOR, NODE)
        self.request_output_variable('R', VECTOR, NODE)
        self.request_output_variable('S', SCALAR, ELEMENT)
        self.request_output_variable('P', SCALAR, ELEMENT)

    def Solve(self):
        """Solves the finite element problem

        """
        # active DOF set dynamically
        self.active_dof = range(self.elements[0].ndof)
        self.validate(ElasticLinknD2, one=True)
        # Assemble the global stiffness and force
        K = self.assemble_global_stiffness()
        F, Q = self.assemble_global_force()
        Kbc, Fbc = self.apply_bc(K, F+Q)
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')

        # Total force, including reaction, and reaction
        Ft = dot(K, self.dofs)
        R = Ft - F - Q

        # reshape R to be the same shape as coord
        u = self.dofs.reshape(self.numnod, -1)
        R = R.reshape(self.numnod, -1)
        p = self.internal_forces(u)
        s = self.stresses(p)
        self.snapshot(U=u, R=R, P=p, S=s)

    # ----------------------------------------------------------------------- #
    # --------------------------- POSTPROCESSING ---------------------------- #
    # ----------------------------------------------------------------------- #
    def internal_forces(self, u):
        """

        .. _truss_int_force:

        Computes the axial internal forces for all elements in the truss

        Parameters
        ----------
        u : array_like
            nodal displacements
            u[i,j] is the jth coordinate displacement of the ith node

        Returns
        -------
        p : ndarray
            Axial internal force

        See Also
        --------
        pyfem2.elemlibN_2.ElasticLinknD2.internal_force

        """
        p = zeros(self.numele)
        for el in self.elements:
            iel = self.mesh.elemap[el.label]
            p[iel] = el.internal_force(u[el.nodes])
        return p

    def stresses(self, p):
        """

        .. _truss_stresses:

        Computes the stress in each truss element

        Parameters
        ----------
        p : ndarray
            Element axial force
            p[e] is the internal force in the eth element

        Returns
        -------
        s : ndarray
            Array of stresses

        """
        return p / array([el.A for el in self.elements])
