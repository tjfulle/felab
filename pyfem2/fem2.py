from numpy import array, dot, zeros
from numpy.linalg import solve, LinAlgError

from .fem1 import FiniteElementModel
from .elemlib import ElasticLink1D2, ElasticLink2D2, ElasticLink3D2
from .constants import *
from .utilities import *

__all__ = ['TrussModel']

class TrussModel(FiniteElementModel):
    dimensions = None
    def Solve(self):
        """Solves the finite element problem

        """
        # ACTIVE DOF SET DYNAMICALLY
        self.active_dof = range(count_digits(self.elements[0].signature))
        self.setup((ElasticLink1D2, ElasticLink2D2, ElasticLink3D2), one=True)

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        du = zeros(self.numdof)
        K, rhs = self.assemble(self.dofs, du)
        Kbc, Fbc = self.apply_bc(K, rhs)
        self.dofs[:] = linsolve(Kbc, Fbc)

        # TOTAL FORCE, INCLUDING REACTION, AND REACTION
        #R = self.assemble(self.dofs, du, cflag=LP_OUTPUT)
        R = dot(K, self.dofs) - rhs

        # RESHAPE R TO BE THE SAME SHAPE AS COORD
        u = self.dofs.reshape(self.numnod, -1)
        R = R.reshape(self.numnod, -1)
        p = self.internal_forces(u)
        s = self.stresses(p)

        frame = self.steps.last.Frame(1.)
        frame.field_outputs['U'].add_data(u)
        frame.field_outputs['P'].add_data(p)
        frame.field_outputs['S'].add_data(s)
        frame.field_outputs['R'].add_data(R)
        frame.converged = True

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
            p[iel] = el.internal_force(u[el.inodes])
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
