import logging
from numpy import *
from numpy.linalg import det, inv

from ..constants import *
from ..utilities import *
from .isop_base import isop_base


class CMDN(isop_base):
    """Continuum stress-displacement element, M dimensions, N nodes"""
    axisymmetric = None
    incompatible_modes = None
    hourglass_control = None
    selective_reduced = None

    @staticmethod
    def variables():
        variables = (('V', SYMTENSOR), ('R', TENSOR), ('E', SYMTENSOR),
                     ('S', SYMTENSOR), ('D', SYMTENSOR))
        variables = (('E', SYMTENSOR), ('DE', SYMTENSOR), ('S', SYMTENSOR),
                     ('V', SYMTENSOR, 1))
        return variables

    def shape(self, *args):
        raise NotImplementedError

    def shapegrad(self, *args):
        raise NotImplementedError

    def response(self, rhs, A, svars, energy, u, du, v, a, time, dtime,
                 kstage, kinc, dltyp, dlmag, predef, lflags,
                 ddlmag, mdload, pnewdt):
        """Assemble the element stiffness and rhs"""

        if lflags[2] not in (STIFF_AND_RHS, STIFF_ONLY, MASS_ONLY, RHS_ONLY, MASS_AND_RHS):
            raise NotImplementedError

        if lflags[2] in (MASS_ONLY, MASS_AND_RHS):
            A[:] = self.mass_matrix(u, du)
            if lflags[2] == MASS_ONLY:
                return

        if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY, MASS_AND_RHS):
            self.eval(rhs, A, svars, u, du, time, dtime, kstage, kinc, predef, lflags)
            if lflags[2] == STIFF_ONLY:
                return

        if kinc == 0:
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
                    logging.warn('UNRECOGNIZED DLTYP FLAG')

    def surface_force(self, edge, qe):

        edgenod = self.edges[edge]
        xb = self.xc[edgenod]

        if self.dimensions == 2:
            if len(xb) == 2:
                # LINEAR SIDE
                gw = ones(2)
                gp = array([-1./sqrt(3.), 1./sqrt(3.)])
                Jac = lambda xi: sqrt((xb[1,1]-xb[0,1])**2+(xb[1,0]-xb[0,0])**2)/2.

            elif len(xb) == 3:
                # QUADRATIC SIDE
                gp = array([-sqrt(3./5.), 0, sqrt(3./5.)])
                gw = array([0.5555555556, 0.8888888889, 0.5555555556])
                def Jac(xi):
                    dxdxi = dot([[-.5 + xi, .5 + xi, -2. * xi]], xb)
                    return sqrt(dxdxi[0, 0] ** 2 + dxdxi[0, 1] ** 2)

            else:
                raise ValueError('UNKNOWN ELEMENT EDGE ORDER')

        else:
            raise ValueError('3D SURFACE FORCE NOT IMPLEMENTED')

        Fe = zeros(self.numdof)
        for (p, xi) in enumerate(gp):
            # FORM GAUSS POINT ON SPECIFIC EDGE
            Ne = self.shape(xi, edge=edge)
            Pe = self.pmatrix(Ne)
            c = Jac(xi) * gw[p]
            if self.axisymmetric == 1:
                rp = dot(Ne, self.xc[:,0])
                c *= rp
            Fe += c * dot(Pe.T, qe)

        return Fe

    def eval(self, rhs, A, svars, u, du, time, dtime, kstage, kinc, predef, lflags):
        """Evaluate the element stiffness and residual"""
        xc = self.xc  # + u.reshape(self.xc.shape)

        if self.incompatible_modes:
            # INCOMPATIBLE MODES STIFFNESSES
            m1 = self.numdofpernod
            Kci = zeros((self.numdof, self.dimensions*m1))
            Kii = zeros((self.dimensions*m1, self.dimensions*m1))

        # DATA FOR INDEXING STATE VARIABLE ARRAY
        ntens = self.ndir + self.nshr
        v = [x[0] for x in self.variables()]
        m = len(v) * ntens
        a1, a2, a3 = [v.index(x) for x in ('E', 'DE', 'S')]

        # COMPUTE INTEGRATION POINT DATA
        for p in range(self.num_gauss):

            # INDEX TO START OF STATE VARIABLES
            ij = m * p
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(xc, xi)

            B = self.bmatrix(dNdx, Ne)

            # STRAIN INCREMENT
            de = dot(B, du)

            # VELOCITY GRADIENT
            # L_ij = dv_i / dx_j = d(du_i/dtime) / dx_j
            #      = du_iI dN_I / dx_j * 1 / dtime
            L = zeros((3, 3))
            L1 = dot(dNdx, du.reshape(-1,self.dimensions))
            L[:self.dimensions,:self.dimensions] = L1 / dtime

            # SYMMETRIC AND DEVIATORIC PARTS -> NEEDED FOR FINITE ROTATIONS
            D = .5 * (L + L.T)
            W = L - D

            V = eye(3)
            I3x3 = eye(3)
            R = eye(3)

            z = -2 * axialv(dot(V, D))
            w = -2. * axialv(W)
            _w_ = w - 2. * dot(inv(V - trace(V) * I3x3), z)
            _W_ = -.5 * axialt(_w_)

            # UPDATE THE ROTATION
            _A = I3x3 - _W_ * dtime / 2.
            _R = dot(I3x3 + _W_ * dtime / 2., R)

            # UPDATED ROTATION
            R = dot(inv(_A), _R)

            # RATE OF STRETCH
            Vdot = dot((D + W), V) - dot(V, _W_)
            V += Vdot * dtime

            # UNROTATE DEFORMATION RATE
            d = dot(R.T, dot(D, R))

            # UNROTATE CAUCHY STRESS
            #T = asmatrix(self.data[0, intpt, STRESS:STRESS+NSYMM])
            #sig = dot(R.T, dot(T, R))

            # CONVERT QUANTITIES TO ARRAYS THAT WILL BE PASSED TO MATERIAL MODEL
            d = asvec(d, self.ndir, self.nshr)
            #de = d*dtime*array([1.,1.,1.,2.])
            #sig = asvec(sig)

            # SET DEFORMATION GRADIENT TO THE IDENTITY
            F0 = eye(self.ndir+self.nshr)
            F = eye(self.ndir+self.nshr)

            # PREDEF AND INCREMENT
            temp = dot(Ne, predef[0,0])
            dtemp = dot(Ne, predef[1,0])

            # MATERIAL RESPONSE
            xv = zeros(1)
            e = array(svars[0,ij+a1*ntens:ij+(a1+1)*ntens])
            s = array(svars[0,ij+a3*ntens:ij+(a3+1)*ntens])
            s, xv, D = self.material.response(
                s, xv, e, de, time, dtime, temp, dtemp, None, None,
                self.ndir, self.nshr, self.ndir+self.nshr, xc, F0, F,
                self.label, kstage, kinc)

            # Rotate stress to material increment
            #T = dot(R, dot(asmatrix(sig), R.T))

            # Calculate strain
            #F = dot(V, R)
            #E = .5 * (dot(F.T, F) - I3x3)

            # STORE THE UPDATED VARIABLES
            svars[1,ij+a1*ntens:ij+(a1+1)*ntens] += de  # STRAIN
            svars[1,ij+a2*ntens:ij+(a2+1)*ntens] = de  # STRAIN INCREMENT
            svars[1,ij+a3*ntens:ij+(a3+1)*ntens] = s  # STRESS

            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            c = J * wt
            if self.axisymmetric == 1:
                rp = dot(Ne, xc[:,0])
                c *= rp

            # Update element residual
            if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
                rhs[:] -=  c * dot(s, B)
                if lflags[2] == RHS_ONLY:
                    continue

            if self.selective_reduced:
                D1, D2 = iso_dev_split(self.ndir, self.nshr, self.dimensions, D)
                A += c * dot(dot(B.T, D2), B)

            else:
                A += c * dot(dot(B.T, D), B)
                if self.axisymmetric == 2:
                    c = J * wt
                    Pe = self.pmatrix(Ne)
                    rp = dot(Ne, xc[:,0])
                    Dh = zeros((2,4))
                    Dh[0,:3] = D[2,:3]
                    Db = zeros((2,4))
                    Db[0,0], Db[1,0] = D[0,0], D[3,0]
                    A += c / rp * dot(Pe.T, (dot(Dh, B) - dot(Db, B)))

            if self.incompatible_modes:
                # INCOMPATIBLE MODES
                G = self.gmatrix(xi)
                Kci += dot(dot(B.T, D), G) * J * wt
                Kii += dot(dot(G.T, D), G) * J * wt

        if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY):
            if self.incompatible_modes:
                A -= dot(dot(Kci, inv(Kii)), Kci.T)

            if self.selective_reduced:
                A += self.sri_correction(svars, u, du, time, dtime,
                                         kstage, kinc, predef, lflags)

            if self.hourglass_control:
                A += self.hourglass_correction(u, du, lflags)

    def hourglass_correction(self, u, du, lflags):

        # PERFORM HOURGLASS CORRECTION
        xc = self.xc
        Khg = zeros((self.numdof, self.numdof))
        for p in range(len(self.hglassp)):

            # SHAPE FUNCTION DERIVATIVE AT HOURGLASS GAUSS POINTS
            xi = array(self.hglassp[p])
            dNdxi = self.shapegrad(xi)

            # JACOBIAN TO NATURAL COORDINATES
            Ne = self.shape(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx, Ne, xi)
            J = det(dxdxi)

            # HOURGLASS BASE VECTORS
            g = array(self.hglassv[p])
            for i in range(len(xi)):
                xi[i] = dot(g, xc[:,i])

            # CORRECT THE BASE VECTORS TO ENSURE ORTHOGONALITY
            scale = 0.
            for a in range(self.nodes):
                for i in range(self.dimensions):
                    g[a] -= xi[i] * dNdx[i,a]
                    scale += dNdx[i,a] * dNdx[i,a]
            scale *= .01 * self.material.G

            for a in range(self.nodes):
                n1 = count_digits(self.signature[a])
                for i in range(n1):
                    for b in range(self.nodes):
                        n2 = count_digits(self.signature[b])
                        for j in range(n2):
                            K = n1 * a + i
                            L = n2 * b + j
                            Khg[K,L] += scale * g[a] * g[b] * J * 4.

        return Khg

    def sri_correction(self, svars, u, du, time, dtime,
                       kstage, kinc, predef, lflags):
        # SELECTIVE REDUCED INTEGRATION CORRECTION

        Ksri = zeros((self.numdof, self.numdof))

        # EVALUATE MATERIAL MODEL AT ELEMENT CENTROID
        xi = self.cp
        xc = self.xc

        # SHAPE FUNCTION AND GRADIENT
        Ne = self.shape(xi)

        # SHAPE FUNCTION DERIVATIVE AT GAUSS POINTS
        dNdxi = self.shapegrad(xi)

        # JACOBIAN TO NATURAL COORDINATES
        dxdxi = dot(dNdxi, xc)
        dxidx = inv(dxdxi)
        J = det(dxdxi)

        # CONVERT SHAPE FUNCTION DERIVATIVES TO DERIVATIVES WRT GLOBAL X
        dNdx = dot(dxidx, dNdxi)
        B = self.bmatrix(dNdx, Ne, xi)

        # STRAIN INCREMENT
        de = dot(B, du)

        # SET DEFORMATION GRADIENT TO THE IDENTITY
        F0 = eye(self.ndir+self.nshr)
        F = eye(self.ndir+self.nshr)

        # PREDEF AND INCREMENT
        temp = dot(Ne, predef[0,0])
        dtemp = dot(Ne, predef[1,0])

        # MATERIAL RESPONSE AT CENTROID
        v = [x[0] for x in self.variables()]
        a1, a2, a3 = [v.index(x) for x in ('E', 'DE', 'S')]
        e = self.interpolate_to_centroid(svars[0], index=a1)
        s = self.interpolate_to_centroid(svars[0], index=a3)
        xv = zeros(1)
        s, xv, D = self.material.response(
            s, xv, e, de, time, dtime, temp, dtemp, None, None,
            self.ndir, self.nshr, self.ndir+self.nshr, xc, F0, F,
            self.label, kstage, kinc)
        D1, D2 = iso_dev_split(self.ndir, self.nshr, self.dimensions, D)

        # GAUSS INTEGRATION
        for p in range(len(self.srip)):
            xi = self.srip[p]
            wt = self.sriw[p]
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            J = det(dxdxi)
            Ne = self.shape(xi)
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx, Ne, xi)
            Ksri += J * wt * dot(dot(B.T, D1), B)

        return Ksri

    def mass_matrix(self, u, du):
        Me = zeros((self.numdof, self.numdof))
        for p in range(self.num_gauss):
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(self.xc, xi)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Pe = self.pmatrix(Ne)
            Me += J * wt * self.material.density * dot(Pe.T, Pe)
        return Me

    def body_force(self, dltyp, dlmag):
        bload = [dlmag[i] for (i, typ) in enumerate(dltyp) if typ==DLOAD]
        Fe = zeros(self.numdof)
        for p in range(self.num_gauss):
            # INDEX TO START OF STATE VARIABLES
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(self.xc, xi)
            Pe = self.pmatrix(Ne)
            for dlmagx in bload:
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                c = J * wt
                if self.axisymmetric == 1:
                    rp = dot(Ne, self.xc[:,0])
                    c *= rp
                Fe += c * dot(Pe.T, dlmagx)
        return Fe
