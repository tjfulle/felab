from numpy import *
from numpy.linalg import inv, det

from .isoplib import (IsoPElement, IsoPReduced, IsoPSelectiveReduced,
                      IsoPIncompatibleModes)

__all__ = ['PlaneStressQuad4', 'PlaneStrainQuad4', 'PlaneStrainQuad4Reduced',
           'PlaneStrainQuad4BBar', 'PlaneStrainQuad4SelectiveReduced',
           'PlaneStressQuad4Incompat']

# --------------------------------------------------------------------------- #
# --------------------- QUADRATIC ISOPARAMETRIC ELEMENTS -------------------- #
# --------------------------------------------------------------------------- #
class IsoPQuad4(IsoPElement):
    """4-node isoparametric element

    Notes
    -----
    Node and element face numbering

              [2]
           3-------2
           |       |
       [3] |       | [1]
           |       |
           0-------1
              [0]

    """
    elefab = {'t':1.}
    signature = (1,1,0,0,0,0,0)
    numdim, numnod, ndof = 2, 4, 2
    gaussp = array([[-1., -1.], [ 1., -1.], [-1.,  1.], [ 1.,  1.]]) / sqrt(3.)
    gaussw = ones(4)
    cp = array([0, 0], dtype=float64)
    xp = array([[-1, -1], [1, -1], [1, 1], [-1, 1],
                [0, -1], [1, 0], [0, 1], [-1, 0]], dtype=float64),
    bgaussw = ones(2)
    bgaussp = array([-1., 1.]) / sqrt(3.)
    edges = array([[0, 1], [1, 2], [2, 3], [3, 0]])

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        A2  = (x[0]*y[1] - x[1]*y[0]) + (x[1]*y[2] - x[2]*y[1])
        A2 += (x[2]*y[3] - x[3]*y[2]) + (x[3]*y[0] - x[0]*y[3])
        return A2 / 2.

    @property
    def volume(self):
        return self.t * self.area

    def isop_map(self, xi):
        raise NotImplementedError

    def shape(self, xi, edge=None):
        if edge is not None:
            # Evaluate shape function on specific edge
            xi = array([[xi,-1.],[1.,xi],[xi,1.],[-1.,xi]][edge])
        N = array([(1. - xi[0]) * (1. - xi[1]),
                   (1. + xi[0]) * (1. - xi[1]),
                   (1. + xi[0]) * (1. + xi[1]),
                   (1. - xi[0]) * (1. + xi[1])]) / 4.
        return N

    def shapegrad(self, xi):
        dN = array([[-1. + xi[1],  1. - xi[1], 1. + xi[1], -1. - xi[1]],
                    [-1. + xi[0], -1. - xi[0], 1. + xi[0],  1. - xi[0]]]) / 4.
        return dN

    def surface_force(self, edge, qe):
        edgenod = self.edges[edge]
        xb = self.xc[edgenod]
        gw = ones(2)
        gp = array([-1./sqrt(3.), 1./sqrt(3.)])
        he = sqrt((xb[1,1]-xb[0,1])**2+(xb[1,0]-xb[0,0])**2)
        Fe = zeros(8)
        for (p, xi) in enumerate(gp):
            # Form Gauss point on specific edge
            Ne = self.shape(xi, edge=edge)
            Pe = self.pmatrix(Ne)
            Fe += he / 2. * gw[p] * dot(Pe.T, qe)
        return Fe

# --------------------------------------------------------------------------- #
# ------------------------ USER ELEMENT TYPES ------------------------------- #
# --------------------------------------------------------------------------- #
class PlaneStressQuad4(IsoPQuad4):
    ndir, nshr = 2, 1
    def bmatrix(self, dN):
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B

class PlaneStrainQuad4(IsoPQuad4):
    ndir, nshr = 3, 1
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B

class PlaneStrainQuad4BBar(IsoPQuad4):
    ndir, nshr = 3, 1
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        # mean dilatational formulation
        dNb = self.shapegradxbar(self.xc)
        for a in range(self.numnod):
            i = 2 * a
            j = i + 1
            bb1 = (dNb[0, a] - dN[0, a]) / 2.
            bb2 = (dNb[1, a] - dN[1, a]) / 2.
            B[0, i:i+2] += [bb1, bb2]
            B[1, i:i+2] += [bb1, bb2]
        return B

class PlaneStrainQuad4Reduced(IsoPQuad4, IsoPReduced):
    ndir, nshr = 3, 1
    gaussp = array([[0., 0.]])
    gaussw = array([4.])
    hglassp = array([[0., 0.]])
    hglassv = array([[1., -1., 1., -1.]]) # hourglass vector
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B

class PlaneStrainQuad4SelectiveReduced(PlaneStrainQuad4, IsoPSelectiveReduced):
    rgaussp = array([[0., 0.]])
    rgaussw = array([4.])

class PlaneStressQuad4Incompat(PlaneStressQuad4, IsoPIncompatibleModes):
    ndir, nshr = 2, 1
    def bmatrix(self, dN):
        """Assemble and return the B matrix"""
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B

    def gmatrix(self, xi):
        """Assemble and return the G matrix"""
        # Algorithm in
        # The Finite Element Method: Its Basis and Fundamentals
        # By Olek C Zienkiewicz, Robert L Taylor, J.Z. Zhu

        # Jacobian at element centroid
        # compute the shape function at the centroid (self.cp)
        dNdxi = self.shapegrad(self.cp)
        
        # compute the deformation gradient at centroid (coordinates are self.xc)
        # and the jacobian
        dxdxi = dot(dNdxi, self.xc)
        dxidx = inv(dxdxi)
        J0 = det(dxidx)
        
        # compute the jacobian of the element
        J = self.jacobian(self.xc, xi)
        
        # compute dNdxi associated with the incompatible modes and then from it
        # and the jacobians computed above compute dNdx
        # N = [1 - xi**2, 1 - eta**2]
        dNdxi = array([[-2*xi[0], 0], [0, -2*xi[1]]])
        dNdx = J0 / J * dot(dxidx, dNdxi)

        # form the incompatible G matrix
        # this is BI in the pdf file
        G = array([[dNdx[0,0], 0,         dNdx[0,1], 0],
                   [0,         dNdx[1,0], 0,         dNdx[1,1]],
                   [dNdx[0,1], dNdx[0,0], dNdx[1,1], dNdx[1,0]]])

        return G
