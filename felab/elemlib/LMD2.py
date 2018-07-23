from numpy import *
from ..utilities import *
from .element_base import element_base


# --------------------------------------------------------------------------- #
# ------------------------------ TRUSS ELEMENT ------------------------------ #
# --------------------------------------------------------------------------- #
class LMD2(element_base):
    """Base class for 2 node elastic link elements

    Parameters
    ----------
    label : int
        Element label
    elenod : list of int
        Internal node IDs of nodes making up this element
    elecoord : ndarray
        Coordinates of nodes
    elemat : object
        A felab.mat.Material instance
    elefab : dict
        Requires area 'A'

    """
    nodes = 2
    elefab = {'A': 1.}

    @classmethod
    def variables(cls):
        return (('P', SCALAR), ('S', SCALAR))

    def response(self, u, du, time, dtime, kstage, kincrement, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, stage_type):
        """Computes the response of a n-dimensional elastic link

        Parameters
        ----------

        Returns
        -------
        k : array_like
            (2*nd, 2*nd) elastic stiffness
        f : array_like
            (2*nd, 0) internal force

        """
        # INTERNAL FORCE
        ndof = count_digits(self.signature[0])

        compute_stiff = cflag in (STIFF_AND_RHS, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_RHS, RHS_ONLY)

        if compute_force:
            fext = zeros(2*ndof)
            fint = zeros(2*ndof)

        if cflag == RHS_ONLY:
            return fext, fint

        # ELEMENT NORMAL
        v = self.xc[1] - self.xc[0]
        h = sqrt(dot(v, v))
        n = v / h
        if self.dimensions == 1:
            nn = 1.
        else:
            nn = outer(n, n)

        svars[1,0] = self.internal_force(u+du)
        svars[1,1] = svars[1,0] / self.A

        if compute_stiff:
            # ASSEMBLE ELEMENT STIFFNESS
            i, j = ndof, 2*ndof
            Ke = zeros((2*ndof, 2*ndof))
            Ke[0:i, 0:i] = Ke[i:j, i:j] =  nn # UPPER LEFT AND LOWER RIGHT 2X2
            Ke[0:i, i:j] = Ke[i:j, 0:i] = -nn # LOWER LEFT AND UPPER RIGHT 2X2
            Ke *= self.A * self.material.E / h

        if cflag == STIFF_AND_RHS:
            return Ke, fext, fint

        elif cflag == STIFF_ONLY:
            return Ke

    def internal_force(self, uc):
        """
        .. _link_int_force:

        Computes the element axial internal force

        Parameters
        ----------
        xc : array_like
            nodal coordinates
            x[i,j] is the jth coordinate of the ith node
        E, A : float
            Young's modulus and cross-sectional area
        uc : array_like
            nodal displacements
            on reshaping to have shape (number of nodes, degrees of freedom),
            u[i,j] is the jth coordinate displacement of the ith node

        Returns
        -------
        p : ndarray
            Array of axial internal forces

        """
        x = self.xc[1] - self.xc[0]
        uc = uc.reshape(self.xc.shape)
        u = uc[1] - uc[0]
        Xu = dot(x, u)
        L = sqrt(dot(x, x))
        return self.material.E * self.A / L * Xu / L
