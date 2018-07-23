from numpy import *
from ..utilities import *
from .DCMDN import DCMDN
from .gauss_rule_info import tri_gauss_rule_info, line_gauss_rule_info


# --------------------------------------------------------------------------- #
# ------------------------ HEAT TRANSFER ELEMENT ---------------------------- #
# --------------------------------------------------------------------------- #
class DC2D3(DCMDN):
    nodes = 3
    dimensions = 2
    num_gauss = 3
    signature = [(0,0,0,0,0,0,1),  # 3 NODE 2D HEAT TRANSFER
                 (0,0,0,0,0,0,1),
                 (0,0,0,0,0,0,1)]
    elefab = {'t': 1.}
    cp = array([1., 1., 1.]) / 3.
    edges = array([[0, 1], [1, 2], [2, 0]])
    xp = array([[1., 0., 0], [0., 1., 0], [0., 0., 1]])

    @staticmethod
    def gauss_rule_info(point):
        return tri_gauss_rule_info(3, point)

    @property
    def area(self):
        x, y = self.xc[:, [0, 1]].T
        a = .5 * (x[0] * (y[1] - y[2]) +
                  x[1] * (y[2] - y[0]) +
                  x[2] * (y[0] - y[1]))
        return a

    @property
    def volume(self):
        return self.t * self.area

    def edge_shape(self, edge, xp):
        # ORDERING OF NODES
        xb = self.xc[self.edges[edge]]
        he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
        o = array({0:[0,1,2],1:[2,0,1],2:[1,2,0]}[edge])
        s = he * (xp + 1) / 2.0
        return array([(he - s) / he, s / he, 0.])[o]

    def shapefun_der(self, coord, qcoord):
        """Shape functions of 3 node triangle

        Parameters
        ----------
        coord : ndarray
            The coordinate in the physical coordinates
        qcoord : ndarray
            The coordinate in the triangle coordinates

        Returns
        -------
        N : ndarray
            The shape function in the natural coordinates
        dN : ndarray
            The shape function derivatives in the physical coordinates
        J : float
            The Jacobian of the transformation

        """
        x = coord[:,0]
        y = coord[:,1]
        z1, z2, z3 = qcoord
        # Triangle coordinates *are* the shape functions
        N = array([z1, z2, z3])
        dNdz = eye(3)
        J = array([[1, 1, 1],
                      [dot(x, dNdz[:,0]), dot(x, dNdz[:,1]), dot(x, dNdz[:,2])],
                      [dot(y, dNdz[:,0]), dot(y, dNdz[:,1]), dot(y, dNdz[:,2])]])
        Jdet = linalg.det(J)
        D = array([[0, 0], [1, 0], [0, 1]])
        P = dot(linalg.inv(J), D)
        dNdx_T = dot(dNdz, P)
        return N, dNdx_T.T, Jdet

    def conduction_stiff_contrib(self):
        # MATERIAL STIFFNESS - "RESISTANCE" TO CONDUCTION
        Ke = zeros((3,3))
        for p in range(self.num_gauss):
            xi, w = self.gauss_rule_info(p)
            Ne, dN, Je  = self.shapefun_der(self.xc, xi)
            k = self.material.model.isotropic_thermal_conductivity(2)
            Ke += Je / 2. * w * dot(dot(dN.T, k), dN)
        return Ke

    def convection_stiff_contrib(self, edge, h):
        # CONVECTION STIFFNESS
        Ke = zeros((3,3))
        for p in range(2):
            xi, w = line_gauss_rule_info(2, p)
            # DETERMINE EDGE LENGTH
            xb = self.xc[self.edges[edge]]
            he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
            s = he * (xi + 1.) / 2.
            N = self.edge_shape(edge, s)
            Ke += h * he / 2. * w * outer(N, N)
        return Ke

    def heat_source(self, f):
        Fe = zeros(3)
        for p in range(self.num_gauss):
            xi, w = self.gauss_rule_info(p)
            Ne, _, Je  = self.shapefun_der(self.xc, xi)
            Fe += Je / 2. * w * Ne * dot(Ne, f)
        return Fe

    def conduction_flux_array(self, edge, qn):
        return self.boundary_flux_array(edge, qn)

    def convection_flux_array(self, edge, Too, h):
        return self.boundary_flux_array(edge, Too * h)

    def boundary_flux_array(self, edge, qn):
        Fe = zeros(3)
        xb = self.xc[self.edges[edge]]
        he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
        for p in range(2):
            xi, w = line_gauss_rule_info(2, p)
            s = he * (xi + 1.) / 2.
            N = self.edge_shape(edge, s)
            Fe += he / 2. * qn * w * N
        return Fe
