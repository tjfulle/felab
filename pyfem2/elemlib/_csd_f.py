import logging
from numpy import *
from numpy.linalg import det, inv

from ._csd import CSDElement
from ..utilities import *

# --------------------------------------------------------------------------- #
# ----------------FULLY INTEGRATED ISOPARAMETRIC ELEMENTS ------------------- #
# --------------------------------------------------------------------------- #
class CSDFElement(CSDElement):
    """Continuum Isoparametric Stress Displacement Element"""
    gaussp = None
    gaussw = None
    variables = ('E', 'DE', 'S')
    integration = None

    def surface_force(self, *args):
        raise NotImplementedError

    def shape(self, *args):
        raise NotImplementedError

    def shapegrad(self, *args):
        raise NotImplementedError

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

        dist = lambda a, b: max(sqrt(dot(a - b, a - b)), 1e-6)
        weights = [1./dist(point, cls.gaussp[i]) for i in range(cls.integration)]

        if data.ndim == 1:
            # SCALAR DATA
            return average(data, weights=weights)

        elif len(data.shape) == 2:
            # VECTOR OR TENSOR DATA
            return average(data, axis=0, weights=weights)

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

            if compute_stiff:
                # ADD CONTRIBUTION OF FUNCTION CALL TO INTEGRAL
                Ke += J * self.gaussw[p] * dot(dot(B.T, D), B)

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
                    # SURFACE LOAD
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
