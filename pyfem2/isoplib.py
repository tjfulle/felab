import logging
from numpy import *
from numpy.linalg import det, inv

# --------------------------------------------------------------------------- #
# ------------------------ ISOPARAMETRIC ELEMENTS --------------------------- #
# --------------------------------------------------------------------------- #
class IsoPElement(object):
    elefab = {}
    ndof, numnod, numdim = None, None, None
    signature = None
    edges = []

    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.nodes = elenod
        self.xc = elecoord
        self.material = elemat
        unknown = [key for key in elefab if key not in self.elefab]
        if unknown:
            raise UserInputError('Unrecognized element fabrication '
                                 'properties: {0}'.format(','.join(unknown)))
        for (name, default) in self.elefab.items():
            p = elefab.get(name, default)
            setattr(self, name, p)

    def surface_force(self, *args):
        raise NotImplementedError

    def pmatrix(self, N):
        S = zeros((self.ndof, self.numnod*self.ndof))
        for i in range(self.numdim):
            S[i, i::self.numdim] = N
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
        dNb = zeros((self.ndof, self.numnod))
        for (npt, xi) in enumerate(self.gaussp):
            # Compute the integrals over the volume
            fac = self.gaussw[npt] * det[npt] / self.numdim / ev
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
        a = zeros((cls.numnod, nx))
        for i in range(cls.numnod):
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
            # Scalar data at each Gauss point of a single element
            assert len(data) == len(cls.gaussp)
            return average(data, weights=weights)

        elif len(data.shape) == 2:
            # Scalar data at each Gauss point
            assert data.shape[1] == len(cls.gaussp)
            return average(data, axis=0, weights=weights)

        elif len(data.shape) == 3:
            assert data.shape[1] == len(cls.gaussp)
            return average(data, axis=1, weights=weights)

        raise TypeError('Unknown data type')

    def stiffness(self, *args):
        """Assemble the element stiffness"""
        # compute integration point data
        Ke = zeros((self.ndof * self.numnod, self.ndof * self.numnod))
        for (p, xi) in enumerate(self.gaussp):
            # Update material state
            J = self.jacobian(self.xc, xi)
            dNdx = self.shapegradx(self.xc, xi)
            B = self.bmatrix(dNdx)
            D = self.material.stiffness(self.ndir, self.nshr)
            # Add contribution of function call to integral
            Ke += dot(dot(B.T, D), B) * J * self.gaussw[p]
        return Ke

    def residual(self, xc, u, stress):
        """Compute the element residual"""
        # compute integration point data
        R = zeros(self.ndof * self.numnod)
        x = self.xc + u
        for (p, xi) in enumerate(self.gaussp):
            J = self.jacobian(xc, xi)
            dNdx = self.shapegradx(xc, xi)
            B = self.bmatrix(dNdx)
            R += dot(stress, B) * J * self.gaussw[p]
        return R

    def force(self, dload, sload):
        """Apply the Von Neumann (natural) BC to the element

        Von Neummann BC is applied as a surface load on element faces.

        Parameters
        ----------
        xc : array_like
            Current nodal coords (of face nodes)

        q : array_like
            Tractions
            trac_vec[i] -> traction on coordinate i as a function of time

        Returns
        -------
        Fe : array_like
            Element distributed load vector

        """
        Fe = zeros(self.ndof * self.numnod)
        for (p, xi) in enumerate(self.gaussp):
            # Body force contribution
            Je = self.jacobian(self.xc, xi)
            Pe = self.pmatrix(self.shape(xi))
            Fe += Je * self.gaussw[p] * dot(Pe.T, dload)

        if not any([dot(qe, qe) >= 1e-12 for qe in sload]):
            return Fe

        for i in range(self.edges.shape[0]):
            # Boundary contribution
            qe = sload[i]
            if dot(qe, qe) <= 1e-12:
                # No force
                continue
            Fe += self.surface_force(i, qe)
        return Fe

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

# --------------------------------------------------------------------------- #
# -------------- REDUCED INTEGRATION ISOPARAMETRIC ELEMENTS ----------------- #
# --------------------------------------------------------------------------- #
class IsoPReduced(IsoPElement):

    def stiffness(self, *args):

        # Get the nominal stiffness
        Kel = super(IsoPReduced, self).stiffness(*args)

        # Perform hourglass correction
        Khg = zeros(Kel.shape)
        for (npt, xi) in enumerate(self.hglassp):
            dN = self.shapegradx(self.xc, xi)
            Je = self.jacobian(self.xc, xi)

            # Hourglass base vectors
            g = self.hglassv[npt]
            for i in range(len(xi)):
                xi[i] = dot(g, self.xc[:,i])

            # Correct the base vectors to ensure orthogonality
            scale = 0.
            for a in range(self.numnod):
                for i in range(self.numdim):
                    g[a] -= xi[i] * dN[i,a]
                    scale += dN[i,a] * dN[i,a]
            scale *= .01 * self.material.G

            for a in range(self.numnod):
                for i in range(self.ndof):
                    for b in range(self.numnod):
                        for j in range(self.ndof):
                            K = self.ndof * a + i
                            L = self.ndof * b + j
                            Khg[K,L] += scale * g[a] * g[b] * Je * 4.
        return Kel + Khg
