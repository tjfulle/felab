import logging
from numpy import *
from numpy.linalg import det, inv

from ...utilities import *
from ..element import Element

class CSDIsoParametricElement(Element):
    """Base class for isoparametric stress-displacement elements"""
    gaussp = None
    gaussw = None
    variables = ('E', 'DE', 'S')
    integration = None
    incompatible_modes = None
    hourglass_control = None
    selective_reduced = None

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

    def pmatrix(self, N):
        n = count_digits(self.signature[0])
        S = zeros((n, self.nodes*n))
        for i in range(self.dimensions):
            S[i, i::self.dimensions] = N
        return S

    @classmethod
    def interpolate_to_centroid(cls, data):
        """Inverse distance weighted average of integration point data at the
        element centroid"""
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

        n = sum([count_digits(nfs) for nfs in self.signature])
        compute_stiff = cflag in (STIFF_AND_RHS, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_RHS, RHS_ONLY, MASS_AND_RHS)
        compute_mass = cflag in (MASS_AND_RHS,)

        if compute_stiff:
            Ke = zeros((n, n))
            if self.incompatible_modes:
                # INCOMPATIBLE MODES STIFFNESSES
                m = count_digits(self.signature[0])
                Kci = zeros((n, self.dimensions*m))
                Kii = zeros((self.dimensions*m, self.dimensions*m))

        if compute_mass:
            Me = zeros((n, n))

        if compute_force:
            xforce = zeros(n)

        if step_type in (GENERAL, DYNAMIC):
            iforce = zeros(n)

        # DATA FOR INDEXING STATE VARIABLE ARRAY
        ntens = self.ndir + self.nshr
        m = len(self.variables) * ntens
        a1, a2, a3 = [self.variables.index(x) for x in ('E', 'DE', 'S')]

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
                    raise NotImplementedError
                else:
                    # STANDARD STIFFNESS
                    Ke += J * self.gaussw[p] * dot(dot(B.T, D), B)

            if compute_mass:
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                Me += J * self.gaussw[p] * self.material.density * outer(N, N)

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
            raise NotImplementedError

        if compute_stiff and self.hourglass_control:

            # PERFORM HOURGLASS CORRECTION
            Khg = zeros(Ke.shape)
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

            Ke += Khg

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

    def surface_force(self, edge, qe):

        edgenod = self.edges[edge]
        xb = self.xc[edgenod]

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

        Fe = zeros(self.numdof)
        for (p, xi) in enumerate(gp):
            # FORM GAUSS POINT ON SPECIFIC EDGE
            Ne = self.shape(xi, edge=edge)
            Pe = self.pmatrix(Ne)
            Fe += Jac(xi) * gw[p] * dot(Pe.T, qe)

        return Fe
