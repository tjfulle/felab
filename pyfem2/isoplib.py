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
    variables = None
    signature = None
    dimensions = None
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
        for (p, xi) in enumerate(self.gaussp):
            # Compute the integrals over the volume
            fac = self.gaussw[p] * det[p] / self.dimensions / ev
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))
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
        ss = 0
        if len(cls.gaussp) == 1:
            weights = [1.]
        else:
            dist = lambda a, b: max(sqrt(dot(a - b, a - b)), 1e-6)
            weights = [1./dist(point, xi) for xi in cls.gaussp]

        if data.ndim == 1:
            # Scalar data
            assert len(data) == len(cls.gaussp)
            return average(data, weights=weights)

        elif len(data.shape) == 2:
            # Vector or tensor data
            assert data.shape[0] == len(cls.gaussp)
            return average(data, axis=0, weights=weights)

        raise TypeError('Unknown data type')

    def update_state(self, u, e, s):
        xc = self.xc  # + u.reshape(self.xc.shape)
        de = zeros_like(e)
        for (p, xi) in enumerate(self.gaussp):
            # UPDATE MATERIAL STATE
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            de[p] = dot(B, u.flatten())
            e[p] += de[p]
            D = self.material.stiffness(self.ndir, self.nshr)
            s[p] += dot(D, de[p])
        return de, e, s

    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload,
                 procedure, nlgeom, cflag, step_type, load_fac):
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

        # COMPUTE INTEGRATION POINT DATA
        bload = [dload[i] for (i, typ) in enumerate(dltyp) if typ==DLOAD]
        for (p, xi) in enumerate(self.gaussp):

            # SHAPE FUNCTION AND GRADIENT
            N = self.shape(xi)

            # SHAPE FUNCTION DERIVATIVE AT GAUSS POINTS
            dNdxi = self.shapegrad(xi)

            # DEFORMATION GRADIENT
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))

            # CONVERT SHAPE FUNCTION DERIVATIVES TO DERIVATIVES WRT GLOBAL X
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)

            # JACOBIAN
            J = det(dxdxi)

            # STRAIN AND INCREMENT
            e = dot(B, u)
            de = dot(B, du)

            # MATERIAL RESPONSE
            s = zeros(self.ndir+self.nshr)
            temp, dtemp = 1., 1.
            xv = zeros(1)
            s, xv, D = self.material.response(s, xv, e, de, time, dtime, temp,
                                              dtemp, self.ndir, self.nshr)

            if compute_stiff:
                # UPDATE MATERIAL STATE
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                Ke += J * self.gaussw[p] * dot(dot(B.T, D), B)

            if compute_force:
                P = self.pmatrix(N)
                for dloadx in bload:
                    # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                    rhs += J * self.gaussw[p] * dot(P.T, dloadx)

            if step_type == GENERAL:
                # SUBTRACT THE RESIDUAL FROM THE FORCE
                resid +=  J * self.gaussw[p] * dot(s, B)

        if cflag == STIFF_ONLY:
            return Ke

        if compute_force:
            for (i, typ) in enumerate(dltyp):
                dloadx = dload[i]
                if typ == DLOAD:
                    continue
                elif typ == SLOAD:
                    # Surface load
                    rhs += self.surface_force(dloadx[0], dloadx[1:])
                else:
                    logging.warn('UNRECOGNIZED DLOAD FLAG')

        rhs *= load_fac

        if step_type == GENERAL:
            rhs -= resid

        if cflag == STIFF_AND_FORCE:
            return Ke, rhs

        elif cflag == FORCE_ONLY:
            return rhs

# --------------------------------------------------------------------------- #
# -------------- REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------------- #
# --------------------------------------------------------------------------- #
class IsoPReduced(IsoPElement):

    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload,
                 procedure, nlgeom, cflag, step_type, load_fac):

        xc = self.xc  # + u.reshape(self.xc.shape)

        # Get the nominal stiffness
        response = super(IsoPReduced, self).response(
            u, du, time, dtime, istep, iframe, dltyp, dload,
            procedure, nlgeom, cflag, step_type, load_fac)

        if cflag == FORCE_ONLY:
            return response

        if cflag == STIFF_AND_FORCE:
            Ke, rhs = response

        elif cflag == STIFF_ONLY:
            Ke = response

        # Perform hourglass correction
        Khg = zeros(Ke.shape)
        for (p, xi) in enumerate(self.hglassp):
            # SHAPE FUNCTION DERIVATIVE AT GAUSS POINTS
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            J = det(dxdxi)

            # Hourglass base vectors
            g = self.hglassv[p]
            for i in range(len(xi)):
                xi[i] = dot(g, xc[:,i])

            # Correct the base vectors to ensure orthogonality
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
# ---------- SELECTIVE REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------- #
# --------------------------------------------------------------------------- #
class IsoPSelectiveReduced(IsoPElement):

    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload,
                 procedure, nlgeom, cflag, step_type, load_fac):
        """Assemble the element stiffness"""

        xc = self.xc  # + u.reshape(self.xc.shape)

        if cflag in (STIFF_AND_FORCE, FORCE_ONLY):
            rhs = super(IsoPSelectiveReduced, self).response(
                u, du, time, dtime, istep, iframe, dltyp, dload,
                procedure, nlgeom, FORCE_ONLY, step_type, load_fac)

            if cflag == FORCE_ONLY:
                return rhs

        # compute integration point data
        n = sum([count_digits(nfs) for nfs in self.signature])
        Ke = zeros((n, n))

        for (p, xi) in enumerate(self.gaussp):
            # SHAPE FUNCTION AND GRADIENT
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            J = det(dxdxi)
            # UPDATE MATERIAL STATE
            D1, D2 = self.material.stiffness(self.ndir, self.nshr, disp=2)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Ke += dot(dot(B.T, D2), B) * J * self.gaussw[p]

        for (p, xi) in enumerate(self.rgaussp):
            # SHAPE FUNCTION AND GRADIENT
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            J = det(dxdxi)
            # UPDATE MATERIAL STATE
            D1, D2 = self.material.stiffness(self.ndir, self.nshr, disp=2)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Ke += dot(dot(B.T, D1), B) * J * self.rgaussw[p]

        if cflag == STIFF_AND_FORCE:
            return Ke, rhs

        elif cflag == STIFF_ONLY:
            return Ke

# --------------------------------------------------------------------------- #
# ---------------- INCOMPATIBLE MODES ISOPARAMETRIC ELEMENTS ---------------- #
# --------------------------------------------------------------------------- #
class IsoPIncompatibleModes(IsoPElement):
    numincomp = None
    def gmatrix(self, xi): raise NotImplementedError
    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload,
                 procedure, nlgeom, cflag, step_type, load_fac):
        """Assemble the element stiffness"""

        xc = self.xc  # + u.reshape(self.xc.shape)

        if cflag in (STIFF_AND_FORCE, FORCE_ONLY):
            rhs = super(IsoPIncompatibleModes, self).response(
                u, du, time, dtime, istep, iframe, dltyp, dload,
                procedure, nlgeom, FORCE_ONLY, step_type, load_fac)

            if cflag == FORCE_ONLY:
                return rhs

        # INCOMPATIBLE MODES STIFFNESSES
        n = sum([count_digits(nfs) for nfs in self.signature])
        m = count_digits(self.signature[0])
        Kcc = zeros((n, n))
        Kci = zeros((n, self.dimensions*m))
        Kii = zeros((self.dimensions*m, self.dimensions*m))

        for (p, xi) in enumerate(self.gaussp):
            # SHAPE FUNCTION GRADIENT
            dNdxi = self.shapegrad(xi)
            dxdxi = dot(dNdxi, xc)
            dxidx = inv(dot(dNdxi, xc))
            dNdx = dot(dxidx, dNdxi)
            B = self.bmatrix(dNdx)
            J = det(dxdxi)
            G = self.gmatrix(xi)
            # UPDATE MATERIAL STATE
            D = self.material.stiffness(self.ndir, self.nshr)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Kcc += dot(dot(B.T, D), B) * J * self.gaussw[p]
            Kci += dot(dot(B.T, D), G) * J * self.gaussw[p]
            Kii += dot(dot(G.T, D), G) * J * self.gaussw[p]
        Ke = Kcc  - dot(dot(Kci, inv(Kii)), Kci.T)

        if cflag == STIFF_AND_FORCE:
            return Ke, rhs

        elif cflag == STIFF_ONLY:
            return Ke
