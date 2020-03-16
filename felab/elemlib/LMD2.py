import numpy as np

from felab.util.numeric import count_digits
from felab.constants import SCALAR, STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY
from felab.elemlib.element_base import element_base


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
    elefab = {"A": 1.0}

    @classmethod
    def variables(cls):
        return (("P", SCALAR), ("S", SCALAR))

    def response(
        self,
        rhs,
        A,
        svars,
        energy,
        u,
        du,
        v,
        a,
        time,
        dtime,
        kstep,
        kinc,
        dltyp,
        dlmag,
        predef,
        lflags,
        ddlmag,
        mdload,
        pnewdt,
    ):
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
        if lflags[2] not in (STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY):
            raise NotImplementedError

        # INTERNAL FORCE
        if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
            rhs[:] = 0.0
            if lflags[2] == RHS_ONLY:
                return

        # ELEMENT NORMAL
        v = self.xc[1] - self.xc[0]
        h = np.sqrt(np.dot(v, v))
        n = v / h
        if self.dimensions == 1:
            nn = 1.0
        else:
            nn = np.outer(n, n)

        svars[1, 0] = self.internal_force(u + du)
        svars[1, 1] = svars[1, 0] / self.A

        if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY):
            ndof = count_digits(self.signature[0])
            i, j = ndof, 2 * ndof
            # ASSEMBLE ELEMENT STIFFNESS
            A[0:i, 0:i] = A[i:j, i:j] = nn  # UPPER LEFT AND LOWER RIGHT 2X2
            A[0:i, i:j] = A[i:j, 0:i] = -nn  # LOWER LEFT AND UPPER RIGHT 2X2
            A *= self.A * self.material.E / h

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
        Xu = np.dot(x, u)
        L = np.sqrt(np.dot(x, x))
        return self.material.E * self.A / L * Xu / L
