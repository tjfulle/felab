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
        self.inodes = elenod
        self.xc = elecoord
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

    def shapegradxbar(self, x):
        det = array([self.jacobian(x, xi) for xi in self.gaussp])
        ev = dot(det, self.gaussw)
        n = sum([count_digits(nfs) for nfs in self.signature])
        dNb = zeros(n).reshape(-1, self.nodes)
        for (npt, xi) in enumerate(self.gaussp):
            # Compute the integrals over the volume
            fac = self.gaussw[npt] * det[npt] / self.dimensions / ev
            dNb += fac*self.shapegradx(x,xi)
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

    def residual(self, xc, stress):
        """Compute the element residual"""
        # compute integration point data
        n = sum([count_digits(nfs) for nfs in self.signature])
        R = zeros(n)
        for (p, xi) in enumerate(self.gaussp):
            J = self.jacobian(xc, xi)
            dNdx = self.shapegradx(xc, xi)
            B = self.bmatrix(dNdx)
            R += dot(stress[p], B) * J * self.gaussw[p]
        return R

    def update_state(self, u, e, s):
        de = zeros_like(e)
        for (p, xi) in enumerate(self.gaussp):
            # Update material state
            J = self.jacobian(self.xc, xi)
            dNdx = self.shapegradx(self.xc, xi)
            B = self.bmatrix(dNdx)
            de[p] = dot(B, u.flatten())
            e[p] += de[p]
            D = self.material.stiffness(self.ndir, self.nshr)
            s[p] += dot(D, de[p])
        return de, e, s

    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload, flags,
                 load_fac):
        """Assemble the element stiffness"""

        x = self.xc + u.reshape(self.xc.shape)

        n = sum([count_digits(nfs) for nfs in self.signature])
        compute_stiff = flags[2] in (1, 2)
        compute_force = flags[2] in (1, 5)
        linear_perturbation = flags[3] == 1

        if compute_stiff:
            Ke = zeros((n, n))

        if compute_force:
            rhs = zeros(n)

        if not linear_perturbation:
            resid = zeros(n)

        # COMPUTE INTEGRATION POINT DATA
        bload = [dload[i] for (i, typ) in enumerate(dltyp) if typ==DLOAD]
        for (p, xi) in enumerate(self.gaussp):

            dNdx = self.shapegradx(self.xc, xi)
            B = self.bmatrix(dNdx)
            J = self.jacobian(self.xc, xi)

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
                P = self.pmatrix(self.shape(xi))
                for dloadx in bload:
                    # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                    rhs += J * self.gaussw[p] * dot(P.T, dloadx)

            if not linear_perturbation:
                # SUBTRACT THE RESIDUAL FROM THE FORCE
                resid +=  J * self.gaussw[p] * dot(s, B)

        if flags[2] == 2:
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

        if not linear_perturbation:
            rhs -= resid

        if flags[2] == 1:
            return Ke, rhs

        elif flags[2] == 5:
            return rhs

# --------------------------------------------------------------------------- #
# -------------- REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------------- #
# --------------------------------------------------------------------------- #
class IsoPReduced(IsoPElement):

    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload, flags,
                 load_fac):

        # Get the nominal stiffness
        response = super(IsoPReduced, self).response(
            u, du, time, dtime, istep, iframe, dltyp, dload, flags, load_fac)
        if flags[2] == 5:
            return response

        if flags[2] == 1:
            Ke, rhs = response

        elif flags[2] == 2:
            Ke = response

        # Perform hourglass correction
        Khg = zeros(Ke.shape)
        for (npt, xi) in enumerate(self.hglassp):
            dN = self.shapegradx(self.xc, xi)
            Je = self.jacobian(self.xc, xi)

            # Hourglass base vectors
            g = self.hglassv[npt]
            for i in range(len(xi)):
                xi[i] = dot(g, self.xc[:,i])

            # Correct the base vectors to ensure orthogonality
            scale = 0.
            for a in range(self.nodes):
                for i in range(self.dimensions):
                    g[a] -= xi[i] * dN[i,a]
                    scale += dN[i,a] * dN[i,a]
            scale *= .01 * self.material.G

            for a in range(self.nodes):
                n1 = count_digits(self.signature[a])
                for i in range(n1):
                    for b in range(self.nodes):
                        n2 = count_digits(self.signature[b])
                        for j in range(n2):
                            K = n1 * a + i
                            L = n2 * b + j
                            Khg[K,L] += scale * g[a] * g[b] * Je * 4.
        Ke += Khg

        if flags[2] == 1:
            return Ke, rhs

        elif flags[2] == 2:
            return Ke

# --------------------------------------------------------------------------- #
# ---------- SELECTIVE REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------- #
# --------------------------------------------------------------------------- #
class IsoPSelectiveReduced(IsoPElement):

    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload, flags,
                 load_fac):
        """Assemble the element stiffness"""

        if flags[2] in (1, 5):
            flags1 = [i for i in flags]
            flags1[2] = 5
            rhs = super(IsoPSelectiveReduced, self).response(
                u, du, time, dtime, istep, iframe, dltyp, dload, flags1, load_fac)

            if flags[2] == 5:
                return rhs

        # compute integration point data
        n = sum([count_digits(nfs) for nfs in self.signature])
        Ke = zeros((n, n))

        for (p, xi) in enumerate(self.gaussp):
            # Update material state
            J = self.jacobian(self.xc, xi)
            dNdx = self.shapegradx(self.xc, xi)
            B = self.bmatrix(dNdx)
            D1, D2 = self.material.stiffness(self.ndir, self.nshr, disp=2)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Ke += dot(dot(B.T, D2), B) * J * self.gaussw[p]

        for (p, xi) in enumerate(self.rgaussp):
            # UPDATE MATERIAL STATE
            J = self.jacobian(self.xc, xi)
            dNdx = self.shapegradx(self.xc, xi)
            B = self.bmatrix(dNdx)
            D1, D2 = self.material.stiffness(self.ndir, self.nshr, disp=2)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Ke += dot(dot(B.T, D1), B) * J * self.rgaussw[p]

        if flags[2] == 1:
            return Ke, rhs

        elif flags[2] == 2:
            return Ke

# --------------------------------------------------------------------------- #
# ---------------- INCOMPATIBLE MODES ISOPARAMETRIC ELEMENTS ---------------- #
# --------------------------------------------------------------------------- #
class IsoPIncompatibleModes(IsoPElement):
    numincomp = None
    def gmatrix(self, xi): raise NotImplementedError
    def response(self, u, du, time, dtime, istep, iframe, dltyp, dload, flags,
                 load_fac):
        """Assemble the element stiffness"""

        if flags[2] in (1, 5):
            flags1 = [i for i in flags]
            flags1[2] = 5
            rhs = super(IsoPIncompatibleModes, self).response(
                u, du, time, dtime, istep, iframe, dltyp, dload, flags1, load_fac)

            if flags[2] == 5:
                return rhs

        # incompatible modes stiffnesses
        n = sum([count_digits(nfs) for nfs in self.signature])
        m = count_digits(self.signature[0])
        Kcc = zeros((n, n))
        Kci = zeros((n, self.dimensions*m))
        Kii = zeros((self.dimensions*m, self.dimensions*m))

        for (p, xi) in enumerate(self.gaussp):
            # UPDATE MATERIAL STATE
            J = self.jacobian(self.xc, xi)
            dNdx = self.shapegradx(self.xc, xi)
            B = self.bmatrix(dNdx)
            G = self.gmatrix(xi)
            D = self.material.stiffness(self.ndir, self.nshr)
            # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
            Kcc += dot(dot(B.T, D), B) * J * self.gaussw[p]
            Kci += dot(dot(B.T, D), G) * J * self.gaussw[p]
            Kii += dot(dot(G.T, D), G) * J * self.gaussw[p]
        Ke = Kcc  - dot(dot(Kci, inv(Kii)), Kci.T)

        if flags[2] == 1:
            return Ke, rhs

        elif flags[2] == 2:
            return Ke
