from numpy import array, dot, zeros
from numpy.linalg import solve, LinAlgError

from .fem1 import FiniteElementModel
from .elemlibN_2 import ElasticLinknD2
from .constants import *

__all__ = ['TrussModel']

class TrussModel(FiniteElementModel):
    numdim = None
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
        e = s / array([el.material.E for el in self.elements])

        frame = self.steps.last.Frame(1.)
        frame.field_outputs['U'].add_data(u)
        frame.field_outputs['E'].add_data(e)
        frame.field_outputs['DE'].add_data(e)
        frame.field_outputs['S'].add_data(s)
        frame.field_outputs['R'].add_data(R)

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
