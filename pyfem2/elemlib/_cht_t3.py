from numpy import *

from ..utilities import *
from ._cht import CHTElement

# --------------------------------------------------------------------------- #
# ------------------------ HEAT TRANSFER ELEMENT ---------------------------- #
# --------------------------------------------------------------------------- #
class PlaneDiffussiveHeatTransferTria3(CHTElement):
    nodes = 3
    signature = [(0,0,0,0,0,0,1),  # 3 NODE 2D HEAT TRANSFER
                 (0,0,0,0,0,0,1),
                 (0,0,0,0,0,0,1)]
    dimensions = 2
    edges = array([[0,1], [1,2], [2,0]])
    gaussp = array([[1., 1.], [4., 1.], [1., 4.]]) / 6.
    gaussw = ones(3) / 3.

    def isop_map(self, xi):
        xc = self.xc
        x = xc[0,0]*xi[0] + xc[1,0]*xi[1] + xc[2,0]*(1-xi[0]-xi[1])
        y = xc[0,1]*xi[0] + xc[1,1]*xi[1] + xc[2,1]*(1-xi[0]-xi[1])
        return x, y

    def jacobian(self):
        # ELEMENT JACOBIAN
        ((X1, Y1), (X2, Y2), (X3, Y3)) = self.xc
        return (X1*Y2 - X2*Y1 - X1*Y3 + X3*Y1 + X2*Y3 - X3*Y2)

    def shape(self, xp):
        Je = self.jacobian()
        ((X1, Y1), (X2, Y2), (X3, Y3)) = self.xc
        x, y = xp
        # ELEMENT SHAPE FUNCTIONS
        Ne = array([X2*Y3 - X3*Y2 + (Y2 - Y3) * x + (X3 - X2) * y,
                    X3*Y1 - X1*Y3 + (Y3 - Y1) * x + (X1 - X3) * y,
                    X1*Y2 - X2*Y1 + (Y1 - Y2) * x + (X2 - X1) * y]) / Je
        return Ne

    def shapegrad(self):
        Je = self.jacobian()
        ((X1, Y1), (X2, Y2), (X3, Y3)) = self.xc
        # ELEMENT SHAPE FUNCTION GRADIENTS
        dN = array([[Y2 - Y3, Y3 - Y1, Y1 - Y2],
                    [X3 - X2, X1 - X3, X2 - X1]]) / Je
        return dN

    def edge_shape(self, edge, xp):
        # ORDERING OF NODES
        xb = self.xc[self.edges[edge]]
        he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
        o = array({0:[0,1,2],1:[2,0,1],2:[1,2,0]}[edge])
        s = he * (xp + 1) / 2.0
        return array([(he - s) / he, s / he, 0.])[o]

    def conduction_stiff_contrib(self):
        # MATERIAL STIFFNESS - "RESISTANCE" TO CONDUCTION
        Je = self.jacobian()
        B = self.shapegrad()
        N = [self.shape(self.isop_map(xi)) for xi in self.gaussp]
        w, k = self.gaussw, self.material.isotropic_thermal_conductivity(2)
        return Je/2.*sum([w[i]*dot(dot(B.T,k),B) for (i, Ni) in enumerate(N)], 0)

    def convection_stiff_contrib(self, edge, h):
        # CONVECTION STIFFNESS
        xi, w = array([-1., 1.])/sqrt(3.), ones(2)
        # DETERMINE EDGE LENGTH
        xb = self.xc[self.edges[edge]]
        he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
        s = he * (xi + 1.) / 2.
        N = [self.edge_shape(edge, si) for si in s]
        return h*he/2.*sum([w[i]*outer(Ni,Ni) for (i, Ni) in enumerate(N)], 0)

    def heat_source(self, f):
        Je = self.jacobian()
        w = self.gaussw
        N = [self.shape(self.isop_map(xi)) for xi in self.gaussp]
        return Je/2.*sum([w[i]*Ni*dot(Ni, f) for (i, Ni) in enumerate(N)],0)

    def conduction_flux_array(self, edge, qn):
        return self.boundary_flux_array(edge, qn)

    def convection_flux_array(self, edge, Too, h):
        return self.boundary_flux_array(edge, Too * h)

    def boundary_flux_array(self, edge, qn):
        xi, w = array([-1., 1.])/sqrt(3.), ones(2)
        xb = self.xc[self.edges[edge]]
        he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
        s = he * (xi + 1.) / 2.
        N = [self.edge_shape(edge, si) for si in s]
        return he/2.*qn*sum([w[i]*Ni for (i, Ni) in enumerate(N)], 0)
