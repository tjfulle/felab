import numpy as np
from numpy.linalg import det, inv

import felab.util.tty as tty
from felab.elemlib.isop_base import isop_base
from felab.util.numeric import asvec, axialt, axialv, count_digits, iso_dev_split
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
    incompatible_modes = None
    hourglass_control = None
    selective_reduced = None

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
            self.eval(rhs, A, svars, u, du, time, dtime, kstep, kframe, predef, lflags)
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

    def eval(self, rhs, A, svars, u, du, time, dtime, kstep, kframe, predef, lflags):
        """Evaluate the element stiffness and residual"""
        xc = self.xc  # + u.reshape(self.xc.shape)

        if self.incompatible_modes:
            # INCOMPATIBLE MODES STIFFNESSES
            m1 = self.numdofpernod
            Kci = np.zeros((self.numdof, self.dimensions * m1))
            Kii = np.zeros((self.dimensions * m1, self.dimensions * m1))

        # DATA FOR INDEXING STATE VARIABLE ARRAY
        ntens = self.ndir + self.nshr
        v = [x[0] for x in self.variables()]
        m = len(v) * ntens
        a1, a2, a3 = [v.index(x) for x in ("E", "DE", "S")]

        # COMPUTE INTEGRATION POINT DATA
        for p in range(self.num_gauss):

            # INDEX TO START OF STATE VARIABLES
            ij = m * p
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(xc, xi)

            B = self.bmatrix(dNdx, Ne)

            # STRAIN INCREMENT
            de = np.dot(B, du)

            # VELOCITY GRADIENT
            # L_ij = dv_i / dx_j = d(du_i/dtime) / dx_j
            #      = du_iI dN_I / dx_j * 1 / dtime
            L = np.zeros((3, 3))
            L1 = np.dot(dNdx, du.reshape(-1, self.dimensions))
            L[: self.dimensions, : self.dimensions] = L1 / dtime

            # SYMMETRIC AND DEVIATORIC PARTS -> NEEDED FOR FINITE ROTATIONS
            D = 0.5 * (L + L.T)
            W = L - D

            V = np.eye(3)
            I3x3 = np.eye(3)
            R = np.eye(3)

            z = -2 * axialv(np.dot(V, D))
            w = -2.0 * axialv(W)
            _w_ = w - 2.0 * np.dot(inv(V - np.trace(V) * I3x3), z)
            _W_ = -0.5 * axialt(_w_)

            # UPDATE THE ROTATION
            _A = I3x3 - _W_ * dtime / 2.0
            _R = np.dot(I3x3 + _W_ * dtime / 2.0, R)

            # UPDATED ROTATION
            R = np.dot(inv(_A), _R)

            # RATE OF STRETCH
            Vdot = np.dot((D + W), V) - np.dot(V, _W_)
            V += Vdot * dtime

            # UNROTATE DEFORMATION RATE
            d = np.dot(R.T, np.dot(D, R))

            # UNROTATE CAUCHY STRESS
            # T = asmatrix(self.data[0, intpt, STRESS:STRESS+NSYMM])
            # sig = np.dot(R.T, np.dot(T, R))

            # CONVERT QUANTITIES TO ARRAYS THAT WILL BE PASSED TO MATERIAL MODEL
            d = asvec(d, self.ndir, self.nshr)
            # de = d*dtime*np.array([1.,1.,1.,2.])
            # sig = asvec(sig)

            # SET DEFORMATION GRADIENT TO THE IDENTITY
            F0 = np.eye(self.ndir + self.nshr)
            F = np.eye(self.ndir + self.nshr)

            # PREDEF AND INCREMENT
            temp = np.dot(Ne, predef[0, 0])
            dtemp = np.dot(Ne, predef[1, 0])

            # MATERIAL RESPONSE
            xv = np.zeros(1)
            e = np.array(svars[0, ij + a1 * ntens : ij + (a1 + 1) * ntens])
            s = np.array(svars[0, ij + a3 * ntens : ij + (a3 + 1) * ntens])
            s, xv, D = self.material.response(
                s,
                xv,
                e,
                de,
                time,
                dtime,
                temp,
                dtemp,
                None,
                None,
                self.ndir,
                self.nshr,
                self.ndir + self.nshr,
                xc,
                F0,
                F,
                self.label,
                kstep,
                kframe,
            )

            # Rotate stress to material increment
            # T = np.dot(R, np.dot(asmatrix(sig), R.T))

            # Calculate strain
            # F = np.dot(V, R)
            # E = .5 * (np.dot(F.T, F) - I3x3)

            # STORE THE UPDATED VARIABLES
            svars[1, ij + a1 * ntens : ij + (a1 + 1) * ntens] += de  # STRAIN
            svars[1, ij + a2 * ntens : ij + (a2 + 1) * ntens] = de  # STRAIN INCREMENT
            svars[1, ij + a3 * ntens : ij + (a3 + 1) * ntens] = s  # STRESS

            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            c = J * wt
            if self.axisymmetric == 1:
                rp = np.dot(Ne, xc[:, 0])
                c *= rp

            # Update element residual
            if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
                rhs[:] -= c * np.dot(s, B)
                if lflags[2] == RHS_ONLY:
                    continue

            if self.selective_reduced:
                D1, D2 = iso_dev_split(self.ndir, self.nshr, self.dimensions, D)
                A += c * np.dot(np.dot(B.T, D2), B)

            else:
                A += c * np.dot(np.dot(B.T, D), B)
                if self.axisymmetric == 2:
                    c = J * wt
                    Pe = self.pmatrix(Ne)
                    rp = np.dot(Ne, xc[:, 0])
                    Dh = np.zeros((2, 4))
                    Dh[0, :3] = D[2, :3]
                    Db = np.zeros((2, 4))
                    Db[0, 0], Db[1, 0] = D[0, 0], D[3, 0]
                    A += c / rp * np.dot(Pe.T, (np.dot(Dh, B) - np.dot(Db, B)))

            if self.incompatible_modes:
                # INCOMPATIBLE MODES
                G = self.gmatrix(xi)
                Kci += np.dot(np.dot(B.T, D), G) * J * wt
                Kii += np.dot(np.dot(G.T, D), G) * J * wt

        if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY):
            if self.incompatible_modes:
                A -= np.dot(np.dot(Kci, inv(Kii)), Kci.T)

            if self.selective_reduced:
                A += self.sri_correction(
                    svars, u, du, time, dtime, kstep, kframe, predef, lflags
                )

            if self.hourglass_control:
                A += self.hourglass_correction(u, du, lflags)

    def hourglass_correction(self, u, du, lflags):

        # PERFORM HOURGLASS CORRECTION
        xc = self.xc
        Khg = np.zeros((self.numdof, self.numdof))
        for p in range(len(self.hglassp)):

            # SHAPE FUNCTION DERIVATIVE AT HOURGLASS GAUSS POINTS
            xi = np.array(self.hglassp[p])
            dNdxi = self.shapegrad(xi)

            # JACOBIAN TO NATURAL COORDINATES
            # Ne = self.shape(xi)
            dxdxi = np.dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            dNdx = np.dot(dxidx, dNdxi)
            # B = self.bmatrix(dNdx, Ne, xi)
            J = det(dxdxi)

            # HOURGLASS BASE VECTORS
            g = np.array(self.hglassv[p])
            for i in range(len(xi)):
                xi[i] = np.dot(g, xc[:, i])

            # CORRECT THE BASE VECTORS TO ENSURE ORTHOGONALITY
            scale = 0.0
            for a in range(self.nodes):
                for i in range(self.dimensions):
                    g[a] -= xi[i] * dNdx[i, a]
                    scale += dNdx[i, a] * dNdx[i, a]
            scale *= 0.01 * self.material.G

            for a in range(self.nodes):
                n1 = count_digits(self.signature[a])
                for i in range(n1):
                    for b in range(self.nodes):
                        n2 = count_digits(self.signature[b])
                        for j in range(n2):
                            K = n1 * a + i
                            L = n2 * b + j
                            Khg[K, L] += scale * g[a] * g[b] * J * 4.0

        return Khg

    def sri_correction(self, svars, u, du, time, dtime, kstep, kframe, predef, lflags):
        # SELECTIVE REDUCED INTEGRATION CORRECTION

        Ksri = np.zeros((self.numdof, self.numdof))

        # EVALUATE MATERIAL MODEL AT ELEMENT CENTROID
        xi = self.cp
        xc = self.xc

        # SHAPE FUNCTION AND GRADIENT
        Ne = self.shape(xi)

        # SHAPE FUNCTION DERIVATIVE AT GAUSS POINTS
        dNdxi = self.shapegrad(xi)

        # JACOBIAN TO NATURAL COORDINATES
        dxdxi = np.dot(dNdxi, xc)
        dxidx = inv(dxdxi)
        J = det(dxdxi)

        # CONVERT SHAPE FUNCTION DERIVATIVES TO DERIVATIVES WRT GLOBAL X
        dNdx = np.dot(dxidx, dNdxi)
        B = self.bmatrix(dNdx, Ne, xi)

        # STRAIN INCREMENT
        de = np.dot(B, du)

        # SET DEFORMATION GRADIENT TO THE IDENTITY
        F0 = np.eye(self.ndir + self.nshr)
        F = np.eye(self.ndir + self.nshr)

        # PREDEF AND INCREMENT
        temp = np.dot(Ne, predef[0, 0])
        dtemp = np.dot(Ne, predef[1, 0])

        # MATERIAL RESPONSE AT CENTROID
        v = [x[0] for x in self.variables()]
        a1, a2, a3 = [v.index(x) for x in ("E", "DE", "S")]
        e = self.interpolate_to_centroid(svars[0], index=a1)
        s = self.interpolate_to_centroid(svars[0], index=a3)
        xv = np.zeros(1)
        s, xv, D = self.material.response(
            s,
            xv,
            e,
            de,
            time,
            dtime,
            temp,
            dtemp,
            None,
            None,
            self.ndir,
            self.nshr,
            self.ndir + self.nshr,
            xc,
            F0,
            F,
            self.label,
            kstep,
            kframe,
        )
        D1, D2 = iso_dev_split(self.ndir, self.nshr, self.dimensions, D)

        # GAUSS INTEGRATION
        for p in range(len(self.srip)):
            xi = self.srip[p]
            wt = self.sriw[p]
            dNdxi = self.shapegrad(xi)
            dxdxi = np.dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            J = det(dxdxi)
            Ne = self.shape(xi)
            dNdx = np.dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx, Ne, xi)
            Ksri += J * wt * np.dot(np.dot(B.T, D1), B)

        return Ksri

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
