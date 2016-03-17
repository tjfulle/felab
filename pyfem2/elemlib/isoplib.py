import logging
from numpy import *
from numpy.linalg import det, inv
from .utilities import *

# --------------------------------------------------------------------------- #
# ------------------------ ISOPARAMETRIC ELEMENTS --------------------------- #
# --------------------------------------------------------------------------- #
class IsoPElement(object):
    elefab = {}
    nodes = None
    variables = ('E', 'DE', 'S')
    signature = None
    dimensions = None
    integration = None
    integration1 = 0
    fully_reduced = False
    selectively_reduced = False
    edges = []

    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.inodes = asarray(elenod)
        self.xc = asarray(elecoord)
        self.material = elemat
        unknown = [key for key in elefab if key not in self.elefab]
        if unknown:
            raise UserInputError('Unrecognized element fabrication '
                                 'properties: {0}'.format(','.join(unknown)))
        for (name, default) in self.elefab.items():
            p = elefab.get(name, default)
            setattr(self, name, p)

        if self.gaussp is None:
            raise NotImplementedError('Gauss points have not been implemented '
                                      'for this element')

    def surface_force(self, *args):
        raise NotImplementedError

    def pmatrix(self, N):
        n = count_digits(self.signature[0])
        S = zeros((n, self.nodes*n))
        for i in range(self.dimensions):
            S[i, i::self.dimensions] = N
        return S

    def jacobian(self, x, xi):
        """Element Jacobian

        1. Shape function derivative at Gauss points
           dNdxi = shapegrad(xi)

        2. Compute the Jacobian matrix J = dNdxi.x
           dxdxi = dot(dNdxi, X)
           J = det(dxdxi)
        """
        dNdxi = self.shapegrad(xi)
        dxdxi = dot(dNdxi, x)
        return det(dxdxi)

    def shapegradx(self, x, xi):
        """Shape function derivatives with respect to global coordinates

        1. shape function derivative at Gauss points
           dNdxi = shapegrad(xi)
           dxdxi = dot(dNdxi, X)

        2. Convert shape function derivatives to derivatives
           wrt global X
           dxidx = inv(dxdxi)
           dNdx = dot(dxidx, dNdxi)
        """
        dNdxi = self.shapegrad(xi)
        dxidx = inv(dot(dNdxi, x))
        return dot(dxidx, dNdxi)

    def shapegradxbar(self, xc):
        det = array([self.jacobian(xc, xi) for xi in self.gaussp])
        ev = dot(det, self.gaussw)
        n = sum([count_digits(nfs) for nfs in self.signature])
        dNb = zeros(n).reshape(-1, self.nodes)
        for p in range(self.integration):
            # COMPUTE THE INTEGRALS OVER THE VOLUME
            xi = self.gaussp[p]
            fac = self.gaussw[p] * det[p] / self.dimensions / ev
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            dNb += fac*dNdx
        return dNb

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

        i1, i2 = cls.integration, 0
        if cls.integration1:
            i1, i2 = cls.integration1, cls.integration-cls.integration1

        if data.ndim == 1:
            # Scalar data
            assert len(data) == i1+i2

        elif len(data.shape) == 2:
            # Vector or tensor data
            assert data.shape[0] == i1+i2

        else:
            raise TypeError('Unknown data type')

        data = data[:i1]

        if i1 == 1:
            weights = [1.]
        else:
            dist = lambda a, b: max(sqrt(dot(a - b, a - b)), 1e-6)
            weights = [1./dist(point, cls.gaussp[i]) for i in range(i1)]

        if data.ndim == 1:
            # Scalar data
            return average(data, weights=weights)

        elif len(data.shape) == 2:
            # Vector or tensor data
            return average(data, axis=0, weights=weights)

        raise TypeError('Unknown data type')

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type, load_fac):
        """Assemble the element stiffness"""

        xc = self.xc  # + u.reshape(self.xc.shape)

        n = sum([count_digits(nfs) for nfs in self.signature])
        compute_stiff = cflag in (STIFF_AND_FORCE, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_FORCE, FORCE_ONLY)

        if compute_stiff:
            Ke = zeros((n, n))

        if compute_force:
            rhs = zeros(n)

        if step_type == GENERAL:
            resid = zeros(n)

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
            s, xv, D = self.material.response(s, xv, e, de, time, dtime, temp,
                                              dtemp, self.ndir, self.nshr, F0, F)

            # STORE THE UPDATED VARIABLES
            svars[1,ij+a1*ntens:ij+(a1+1)*ntens] += de  # STRAIN
            svars[1,ij+a2*ntens:ij+(a2+1)*ntens] = de  # STRAIN INCREMENT
            svars[1,ij+a3*ntens:ij+(a3+1)*ntens] = s  # STRESS

            if self.selectively_reduced:
                # SELECTIVELY REDUCED INTEGRATION
                D1, D2 = iso_dev_split(self.ndir, self.nshr, D)
                D = D1 if p >= self.integration1 else D2

            if compute_stiff:
                if self.fully_reduced and p < self.integration1:
                    pass
                else:
                    # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                    Ke += J * self.gaussw[p] * dot(dot(B.T, D), B)

            if self.integration1 and p >= self.integration1:
                # DO NOT COMPUTE FORCE FOR REDUCED INTEGRATED PORTION
                continue

            if compute_force:
                P = self.pmatrix(N)
                for dloadx in bload:
                    # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                    rhs += J * self.gaussw[p] * dot(P.T, dloadx)

            if step_type == GENERAL:
                # UPDATE THE RESIDUAL
                resid +=  J * self.gaussw[p] * dot(s, B)

        if cflag == LP_OUTPUT:
            return

        if cflag == STIFF_ONLY:
            return Ke

        if compute_force:
            for (i, typ) in enumerate(dltyp):
                if typ == DLOAD:
                    continue
                dloadx = dload[i]
                if typ == SLOAD:
                    # Surface load
                    rhs += self.surface_force(dloadx[0], dloadx[1:])
                else:
                    logging.warn('UNRECOGNIZED DLOAD FLAG')

        if compute_force:
            rhs *= load_fac

        if step_type == GENERAL:
            # SUBTRACT RESIDUAL FROM INTERNAL FORCE
            rhs -= resid

        if cflag == STIFF_AND_FORCE:
            return Ke, rhs

        elif cflag == FORCE_ONLY:
            return rhs

# --------------------------------------------------------------------------- #
# -------------- REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------------- #
# --------------------------------------------------------------------------- #
class IsoPReduced(IsoPElement):

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type, load_fac):

        xc = self.xc  # + u.reshape(self.xc.shape)

        # GET THE NOMINAL RESPONSE
        response = super(IsoPReduced, self).response(
            u, du, time, dtime, istep, iframe, svars, dltyp, dload,
            predef, procedure, nlgeom, cflag, step_type, load_fac)

        if cflag == LP_OUTPUT:
            return

        elif cflag == FORCE_ONLY:
            return response

        elif cflag == STIFF_AND_FORCE:
            Ke, rhs = response

        elif cflag == STIFF_ONLY:
            Ke = response

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

        if cflag == STIFF_AND_FORCE:
            return Ke, rhs

        elif cflag == STIFF_ONLY:
            return Ke

# --------------------------------------------------------------------------- #
# ---------------- INCOMPATIBLE MODES ISOPARAMETRIC ELEMENTS ---------------- #
# --------------------------------------------------------------------------- #
class IsoPIncompatibleModes(IsoPElement):

    def gmatrix(self, xi):
        raise NotImplementedError

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type, load_fac):
        """Assemble the element stiffness"""

        xc = self.xc  # + u.reshape(self.xc.shape)

        response = super(IsoPIncompatibleModes, self).response(
            u, du, time, dtime, istep, iframe, svars, dltyp, dload,
            predef, procedure, nlgeom, cflag, step_type, load_fac)

        if cflag == LP_OUTPUT:
            return

        elif cflag == FORCE_ONLY:
            return response

        elif cflag == STIFF_AND_FORCE:
            Kcc, rhs = response

        elif cflag == STIFF_ONLY:
            Kcc = response

        # INCOMPATIBLE MODES STIFFNESSES
        n = Kcc.shape[0]
        m = count_digits(self.signature[0])
        Kci = zeros((n, self.dimensions*m))
        Kii = zeros((self.dimensions*m, self.dimensions*m))

        # DATA FOR INDEXING STATE VARIABLE ARRAY
        ntens = self.ndir + self.nshr
        m = len(self.variables) * ntens
        a1, a2, a3 = [self.variables.index(x) for x in ('E', 'DE', 'S')]

        # COMPUTE INTEGRATION POINT DATA
        for (p, xi) in enumerate(self.gaussp):
            ij = m * p

            # SHAPE FUNCTION GRADIENT
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dxdxi)
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            J = det(dxdxi)
            G = self.gmatrix(xi)

            # STRAIN AND INCREMENT
            de = dot(B, du)

            # SET DEFORMATION GRADIENT TO THE IDENTITY
            F0 = eye(self.ndir+self.nshr)
            F = eye(self.ndir+self.nshr)

            # PREDEF AND INCREMENT
            N = self.shape(xi)
            temp = dot(N, predef[0,0])
            dtemp = dot(N, predef[1,0])

            # MATERIAL RESPONSE
            xv = zeros(1)
            e = svars[0,ij+a1*ntens:ij+(a1+1)*ntens]
            s = svars[0,ij+a3*ntens:ij+(a3+1)*ntens]
            s, xv, D = self.material.response(s, xv, e, de, time, dtime, temp,
                                              dtemp, self.ndir, self.nshr, F0, F)

            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Kci += dot(dot(B.T, D), G) * J * self.gaussw[p]
            Kii += dot(dot(G.T, D), G) * J * self.gaussw[p]

        Ke = Kcc  - dot(dot(Kci, inv(Kii)), Kci.T)

        if cflag == STIFF_AND_FORCE:
            return Ke, rhs

        elif cflag == STIFF_ONLY:
            return Ke
