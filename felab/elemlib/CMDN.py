import numpy as np

import felab.util.tty as tty
from felab.elemlib.isop_base import isop_base
from felab.constants import (
    SYMTENSOR,
    TENSOR,
    STIFF_AND_RHS,
    STIFF_ONLY,
    MASS_ONLY,
    MASS_AND_RHS,
    DLOAD,
    SLOAD,
    RHS_ONLY,
)


class CMDN(isop_base):
    """Continuum stress-displacement element, M dimensions, N nodes"""
    axisymmetric = None

    @staticmethod
    def variables():
        variables = (
            ("V", SYMTENSOR),
            ("R", TENSOR),
            ("E", SYMTENSOR),
            ("S", SYMTENSOR),
            ("D", SYMTENSOR),
        )
        variables = (
            ("E", SYMTENSOR),
            ("DE", SYMTENSOR),
            ("S", SYMTENSOR),
            ("V", SYMTENSOR, 1),
        )
        return variables

    def shape(self, *args):
        raise NotImplementedError

    def shapegrad(self, *args):
        raise NotImplementedError

    def eval(
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
        kframe,
        dltyp,
        dlmag,
        predef,
        lflags,
        ddlmag,
        mdload,
        pnewdt,
    ):
        """Assemble the element stiffness and rhs"""

        if lflags[2] not in (
            STIFF_AND_RHS,
            STIFF_ONLY,
            MASS_ONLY,
            RHS_ONLY,
            MASS_AND_RHS,
        ):
            raise NotImplementedError

        if lflags[2] in (MASS_ONLY, MASS_AND_RHS):
            A[:] = self.mass_matrix(u, du)
            if lflags[2] == MASS_ONLY:
                return

        if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY, MASS_AND_RHS):
            A += self.stiffness(svars, u, du, time, dtime, kstep, kframe, predef)
            if lflags[2] == STIFF_ONLY:
                return

        if kframe == 0:
            return

        if lflags[2] in (STIFF_AND_RHS, RHS_ONLY, MASS_AND_RHS):
            rhs[:] += self.body_force(dltyp, dlmag)
            for (i, typ) in enumerate(dltyp):
                if typ == DLOAD:
                    continue
                if typ == SLOAD:
                    # SURFACE LOAD
                    iedge, components = dlmag[i][0], dlmag[i][1:]
                    rhs[:] += self.surface_force(iedge, components)
                else:
                    tty.warn("UNRECOGNIZED DLTYP FLAG")

    def surface_force(self, edge, qe):

        edgenod = self.edges[edge]
        xb = self.xc[edgenod]

        if self.dimensions == 2:
            if len(xb) == 2:
                # LINEAR SIDE
                gw = np.ones(2)
                gp = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
                Jac = (
                    lambda xi: np.sqrt(
                        (xb[1, 1] - xb[0, 1]) ** 2 + (xb[1, 0] - xb[0, 0]) ** 2
                    )
                    / 2.0
                )

            elif len(xb) == 3:
                # QUADRATIC SIDE
                gp = np.array([-np.sqrt(3.0 / 5.0), 0, np.sqrt(3.0 / 5.0)])
                gw = np.array([0.5555555556, 0.8888888889, 0.5555555556])

                def Jac(xi):
                    dxdxi = np.dot([[-0.5 + xi, 0.5 + xi, -2.0 * xi]], xb)
                    return np.sqrt(dxdxi[0, 0] ** 2 + dxdxi[0, 1] ** 2)

            else:
                raise ValueError("UNKNOWN ELEMENT EDGE ORDER")

        else:
            raise ValueError("3D SURFACE FORCE NOT IMPLEMENTED")

        Fe = np.zeros(self.numdof)
        for (p, xi) in enumerate(gp):
            # FORM GAUSS POINT ON SPECIFIC EDGE
            Ne = self.shape(xi, edge=edge)
            Pe = self.pmatrix(Ne)
            c = Jac(xi) * gw[p]
            if self.axisymmetric == 1:
                rp = np.dot(Ne, self.xc[:, 0])
                c *= rp
            Fe += c * np.dot(Pe.T, qe)

        return Fe

    def mass_matrix(self, u, du):
        Me = np.zeros((self.numdof, self.numdof))
        for p in range(self.num_gauss):
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(self.xc, xi)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Pe = self.pmatrix(Ne)
            Me += J * wt * self.material.density * np.dot(Pe.T, Pe)
        return Me

    def body_force(self, dltyp, dlmag):
        bload = [dlmag[i] for (i, typ) in enumerate(dltyp) if typ == DLOAD]
        Fe = np.zeros(self.numdof)
        for p in range(self.num_gauss):
            # INDEX TO START OF STATE VARIABLES
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(self.xc, xi)
            Pe = self.pmatrix(Ne)
            for dlmagx in bload:
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                c = J * wt
                if self.axisymmetric == 1:
                    rp = np.dot(Ne, self.xc[:, 0])
                    c *= rp
                Fe += c * np.dot(Pe.T, dlmagx)
        return Fe
