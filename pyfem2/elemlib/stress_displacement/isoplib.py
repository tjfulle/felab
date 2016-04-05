import logging
from numpy import *
from numpy.linalg import det, inv

from ...constants import *
from ...utilities import *
from ..element import Element

class CSDIsoParametricElement(Element):
    """Base class for isoparametric stress-displacement elements"""
    gaussp = None
    gaussw = None
    integration = None
    incompatible_modes = None
    hourglass_control = None
    selective_reduced = None

    @classmethod
    def variables(self):
        variables = (('V', SYMTENSOR), ('R', TENSOR), ('E', SYMTENSOR),
                     ('S', SYMTENSOR), ('D', SYMTENSOR))
        variables = (('E', SYMTENSOR), ('DE', SYMTENSOR), ('S', SYMTENSOR),
                     ('V', SYMTENSOR, 1))
        return variables

    def shape(self, *args):
        raise NotImplementedError

    def shapegrad(self, *args):
        raise NotImplementedError

    def gmatrix(self, *args):
        return NotImplementedError

    def bmatrix(self, *args):
        return NotImplementedError

    @property
    def numdof(self):
        return sum([count_digits(nfs) for nfs in self.signature])
    @property
    def numdofpernod(self):
        return count_digits(self.signature[0])

    def pmatrix(self, N):
        n = count_digits(self.signature[0])
        S = zeros((n, self.nodes*n))
        for i in range(self.dimensions):
            S[i, i::self.dimensions] = N
        return S

    @classmethod
    def interpolate_to_centroid(cls, data, index=None):
        """Inverse distance weighted average of integration point data at the
        element centroid"""
        if index is not None:
            ntens = cls.ndir + cls.nshr
            m = len(cls.variables()) * ntens
            data = row_stack([data[(m*p)+index*ntens:(m*p)+(index+1)*ntens]
                              for p in range(cls.integration)])
        return cls.average(cls.cp, data)

    @classmethod
    def project_to_nodes(cls, data, v):
        """Inverse distance weighted average of integration point data at each
        element node"""
        nx = len(v)
        a = zeros((cls.nodes, nx))
        for i in range(cls.nodes):
            a[i,:] = cls.average(cls.xp[i], data, v)
        return a

    @classmethod
    def average(cls, point, data, v=None):
        """Inverse distance weighted average of integration point data at point"""

        if data.ndim == 1:
            # SCALAR DATA
            assert len(data) == cls.integration

        elif len(data.shape) == 2:
            # VECTOR OR TENSOR DATA
            assert data.shape[0] == cls.integration

        else:
            raise TypeError('Unknown data type')

        if cls.integration == 1:
            weights = [1.]
        else:
            dist = lambda a, b: max(sqrt(dot(a - b, a - b)), 1e-6)
            weights = [1./dist(point,cls.gaussp[i]) for i in range(cls.integration)]

        if data.ndim == 1:
            # SCALAR DATA
            return average(data, weights=weights)

        elif len(data.shape) == 2:
            # VECTOR OR TENSOR DATA
            return average(data, axis=0, weights=weights)

    def response(self, u, du, time, dtime, kstep, kframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type):
        """Assemble the element stiffness and rhs"""

        xc = self.xc  # + u.reshape(self.xc.shape)

        compute_stiff = cflag in (STIFF_AND_RHS, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_RHS, RHS_ONLY, MASS_AND_RHS)
        compute_mass = cflag in (MASS_AND_RHS,)

        n = self.numdof
        if compute_stiff:
            Ke = zeros((n, n))
            if self.incompatible_modes:
                # INCOMPATIBLE MODES STIFFNESSES
                m1 = self.numdofpernod
                Kci = zeros((n, self.dimensions*m1))
                Kii = zeros((self.dimensions*m1, self.dimensions*m1))

        if compute_mass:
            Me = zeros((n, n))

        if compute_force:
            xforce = zeros(n)

        if step_type in (GENERAL, DYNAMIC):
            iforce = zeros(n)

        # DATA FOR INDEXING STATE VARIABLE ARRAY
        ntens = self.ndir + self.nshr
        v = [x[0] for x in self.variables()]
        m = len(v) * ntens
        a1, a2, a3 = [v.index(x) for x in ('E', 'DE', 'S')]

        # COMPUTE INTEGRATION POINT DATA
        bload = [dload[i] for (i, typ) in enumerate(dltyp) if typ==DLOAD]
        for p in range(self.integration):

            # INDEX TO START OF STATE VARIABLES
            ij = m * p

            # SHAPE FUNCTION AND GRADIENT
            xi = self.gaussp[p]
            N = self.shape(xi)

            # SHAPE FUNCTION DERIVATIVE AT GAUSS POINTS
            dNdxi = self.shapegrad(xi)

            # JACOBIAN TO NATURAL COORDINATES
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            J = det(dxdxi)

            # CONVERT SHAPE FUNCTION DERIVATIVES TO DERIVATIVES WRT GLOBAL X
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)

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
            A = I3x3 - _W_ * dtime / 2.
            RHS = dot(I3x3 + _W_ * dtime / 2., R)

            # UPDATED ROTATION
            R = dot(inv(A), RHS)

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
            temp = dot(N, predef[0,0])
            dtemp = dot(N, predef[1,0])

            # MATERIAL RESPONSE
            xv = zeros(1)
            e = svars[0,ij+a1*ntens:ij+(a1+1)*ntens]
            s = svars[0,ij+a3*ntens:ij+(a3+1)*ntens]
            s, xv, D = self.material.response(
                s, xv, e, de, time, dtime, temp, dtemp, None, None,
                self.ndir, self.nshr, self.ndir+self.nshr, xc, F0, F,
                self.label, kstep, kframe)

            # Rotate stress to material frame
            #T = dot(R, dot(asmatrix(sig), R.T))

            # Calculate strain
            #F = dot(V, R)
            #E = .5 * (dot(F.T, F) - I3x3)

            # STORE THE UPDATED VARIABLES
            svars[1,ij+a1*ntens:ij+(a1+1)*ntens] += de  # STRAIN
            svars[1,ij+a2*ntens:ij+(a2+1)*ntens] = de  # STRAIN INCREMENT
            svars[1,ij+a3*ntens:ij+(a3+1)*ntens] = s  # STRESS

            if compute_stiff:
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                if self.incompatible_modes:
                    # INCOMPATIBLE MODES
                    G = self.gmatrix(xi)
                    Kci += dot(dot(B.T, D), G) * J * self.gaussw[p]
                    Kii += dot(dot(G.T, D), G) * J * self.gaussw[p]
                elif self.selective_reduced:
                    D1, D2 = iso_dev_split(self.ndir, self.nshr, self.dimensions, D)
                    Ke += J * self.gaussw[p] * dot(dot(B.T, D2), B)
                else:
                    # STANDARD STIFFNESS
                    Ke += J * self.gaussw[p] * dot(dot(B.T, D), B)

            if compute_mass:
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                P = self.pmatrix(N)
                Me += J * self.gaussw[p] * self.material.density * dot(P.T, P)

            if compute_force:
                P = self.pmatrix(N)
                for dloadx in bload:
                    # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                    xforce += J * self.gaussw[p] * dot(P.T, dloadx)

            if step_type == GENERAL:
                # UPDATE THE RESIDUAL
                iforce +=  J * self.gaussw[p] * dot(s, B)

        if cflag == LP_OUTPUT:
            return

        if compute_stiff and self.incompatible_modes:
            Ke -= dot(dot(Kci, inv(Kii)), Kci.T)

        if compute_stiff and self.selective_reduced:
            Ke += self.sri_correction(u, du, time, dtime, kstep, kframe,
                                      svars, predef, nlgeom)

        if compute_stiff and self.hourglass_control:
            Ke += self.hourglass_correction(u, du, nlgeom)

        if cflag == STIFF_ONLY:
            return Ke

        if cflag == MASS_ONLY:
            return Me

        if compute_force:
            for (i, typ) in enumerate(dltyp):
                if typ == DLOAD:
                    continue
                if typ == SLOAD:
                    # SURFACE LOAD
                    iedge, components = dload[i][0], dload[i][1:]
                    xforce += self.surface_force(iedge, components)
                else:
                    logging.warn('UNRECOGNIZED DLOAD FLAG')

        if step_type in (GENERAL, DYNAMIC) and compute_force:
            # SUBTRACT RESIDUAL FROM INTERNAL FORCE
            rhs = xforce - iforce

        else:
            rhs = xforce

        if cflag == STIFF_AND_RHS:
            return Ke, rhs

        elif cflag == RHS_ONLY:
            return rhs

        elif cflag == MASS_AND_RHS:
            return Me, rhs

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
            Fe += Jac(xi) * gw[p] * dot(Pe.T, qe)

        return Fe

    def hourglass_correction(self, u, du, nlgeom):

        # PERFORM HOURGLASS CORRECTION
        xc = self.xc
        n = self.numdof
        Khg = zeros((n, n))
        for p in range(len(self.hglassp)):

            # SHAPE FUNCTION DERIVATIVE AT HOURGLASS GAUSS POINTS
            xi = array(self.hglassp[p])
            dNdxi = self.shapegrad(xi)

            # JACOBIAN TO NATURAL COORDINATES
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
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

    def sri_correction(self, u, du, time, dtime, kstep, kframe, svars, predef,
                       nlgeom):
        # SELECTIVE REDUCED INTEGRATION CORRECTION

        n = self.numdof
        Ksri = zeros((n, n))

        # EVALUATE MATERIAL MODEL AT ELEMENT CENTROID
        xi = self.cp
        xc = self.xc

        # SHAPE FUNCTION AND GRADIENT
        N = self.shape(xi)

        # SHAPE FUNCTION DERIVATIVE AT GAUSS POINTS
        dNdxi = self.shapegrad(xi)

        # JACOBIAN TO NATURAL COORDINATES
        dxdxi = dot(dNdxi, xc)
        dxidx = inv(dxdxi)
        J = det(dxdxi)

        # CONVERT SHAPE FUNCTION DERIVATIVES TO DERIVATIVES WRT GLOBAL X
        dNdx = dot(dxidx, dNdxi)
        B = self.bmatrix(dNdx)

        # STRAIN INCREMENT
        de = dot(B, du)

        # SET DEFORMATION GRADIENT TO THE IDENTITY
        F0 = eye(self.ndir+self.nshr)
        F = eye(self.ndir+self.nshr)

        # PREDEF AND INCREMENT
        temp = dot(N, predef[0,0])
        dtemp = dot(N, predef[1,0])

        # MATERIAL RESPONSE AT CENTROID
        v = [x[0] for x in self.variables()]
        a1, a2, a3 = [v.index(x) for x in ('E', 'DE', 'S')]
        e = self.interpolate_to_centroid(svars[0], index=a1)
        s = self.interpolate_to_centroid(svars[0], index=a3)
        xv = zeros(1)
        s, xv, D = self.material.response(
            s, xv, e, de, time, dtime, temp, dtemp, None, None,
            self.ndir, self.nshr, self.ndir+self.nshr, xc, F0, F,
            self.label, kstep, kframe)
        D1, D2 = iso_dev_split(self.ndir, self.nshr, self.dimensions, D)

        # GAUSS INTEGRATION
        for p in range(len(self.srip)):
            xi = self.srip[p]
            w = self.sriw[p]
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            J = det(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            Ksri += J * w * dot(dot(B.T, D1), B)

        return Ksri
