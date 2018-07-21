from numpy import *
from ..utilities import *
from .isop_base import heat_transfer
from .isop_p3_base import isop_p3_base
from .gauss_rule_info import tri_gauss_rule_info, line_gauss_rule_info


# --------------------------------------------------------------------------- #
# ------------------------ HEAT TRANSFER ELEMENT ---------------------------- #
# --------------------------------------------------------------------------- #
class DC2D3(isop_p3_base, heat_transfer):
    num_gauss = 3
    signature = [(0,0,0,0,0,0,1),  # 3 NODE 2D HEAT TRANSFER
                 (0,0,0,0,0,0,1),
                 (0,0,0,0,0,0,1)]

    @staticmethod
    def gauss_rule_info(point):
        return tri_gauss_rule_info(3, point)

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
