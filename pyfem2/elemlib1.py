from numpy import *
from utilities import *

__all__ = ['LinknD2', 'Link1D2', 'Link2D2', 'Link3D2',
           'BeamColumn2D', 'DiffussiveHeatTransfer2D3', 'Tria3', 'Quad4']

# --------------------------------------------------------------------------- #
# --------------------------- ELEMENT LIBRARY ------------------------------- #
# --------------------------------------------------------------------------- #
def ElementFamily(numdim, numnod):
    if numnod == 2:
        return LinknD2
    elif numdim == 2 and numnod == 3:
        return Tria3
    elif numdim == 2 and numnod == 6:
        return Tria6
    elif numdim == 2 and numnod == 4:
        return Quad4
    elif numdim == 2 and numnod == 8:
        return Quad8
    raise ValueError('Unknown element family')

# --------------------------------------------------------------------------- #
# -------------------------- BASE ELMENT CLASS ------------------------------ #
# --------------------------------------------------------------------------- #
class Element(object):
    nfab = 0
    ndof, numnod, numdim = None, None, None
    gaussp, gaussw = None, None
    signature = None
    edges = []
    def stiffness(self, *args): raise NotImplementedError
    def force(self, *args): raise NotImplementedError

# --------------------------------------------------------------------------- #
# ------------------------------ TRUSS ELEMENT ------------------------------ #
# --------------------------------------------------------------------------- #
class LinknD2(Element):
    nfab = 1
    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.nodes = asarray(elenod, dtype=int)
        self.xc = asarray(elecoord, dtype=float)
        self.material = elemat
        self.A = elefab.get('A')
        if self.A is None:
            raise ValueError('Incorrect number of element fabrication properties')

    def stiffness(self, *args):
        """Computes the stiffness of a n-dimensional elastic link

        Parameters
        ----------
        xc : array_like
            nodal coordinates
            xc[i,j] is the jth coordinate of the ith node
        E, A : float
            Young's modulus and area of the bar

        Returns
        -------
        k : array_like
            2*nd x 2*nd elastic stiffness

        """
        # Element normal
        v = self.xc[1] - self.xc[0]
        h = sqrt(dot(v, v))
        n = v / h
        if self.numdim == 1:
            nn = 1.
        else:
            nn = outer(n, n)

        # Assemble element stiffness
        i, j = self.ndof, 2*self.ndof
        k = zeros((2*self.ndof, 2*self.ndof))
        k[0:i, 0:i] = k[i:j, i:j] =  nn # upper left and lower right 2x2
        k[0:i, i:j] = k[i:j, 0:i] = -nn # lower left and upper right 2x2
        return self.A * self.material.E / h * k

    def force(self, *args):
        return zeros(self.numnod*self.ndof)

    def internal_force(self, uc):
        """Computes the element axial internal force

        Parameters
        ----------
        xc : array_like
            nodal coordinates
            x[i,j] is the jth coordinate of the ith node
        E, A : float
            Young's modulus and cross-sectional area
        uc : array_like
            nodal displacements
            on reshaping to have shape (number of nodes, degrees of freedom),
            u[i,j] is the jth coordinate displacement of the ith node

        Returns
        -------
        p : ndarray
            Array of axial internal forces

        """
        x = self.xc[1] - self.xc[0]
        u = uc[1] - uc[0]
        Xu = dot(x, u)
        L = sqrt(dot(x, x))
        return self.material.E * self.A / L * Xu / L

class Link1D2(LinknD2):
    numdim, ndof, numnod = 1, 1, 2
    signature = 1000000  # 2 NODE 1D LINE LINK

class Link2D2(LinknD2):
    numdim, ndof, numnod = 2, 2, 2
    signature = 1100000  # 2 NODE 2D LINE LINK

class Link3D2(LinknD2):
    numdim, ndof, numnod = 3, 3, 2
    signature = 1110000  # 2 NODE 3D LINE LINK

# --------------------------------------------------------------------------- #
# ---------------------- EULER BERNOULI BEAM ELEMENT ------------------------ #
# --------------------------------------------------------------------------- #
class BeamColumn2D(Element):
    nfab = 2
    numdim, ndof, numnod = 2, 3, 2
    signature = 1100010

    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.nodes = asarray(elenod, dtype=int)
        self.xc = asarray(elecoord, dtype=float)
        self.material = elemat
        self.A, self.Izz = elefab.get('A'), elefab.get('Izz')
        if self.A is None or self.Izz is None:
            raise ValueError('Incorrect element fabrication properties')

    def stiffness(self, *args):
        # Compute element normal
        v = self.xc[1] - self.xc[0]
        h = sqrt(dot(v, v))
        n = v / h
        # Transformation matrix
        Te = eye(6)
        Te[0:2, 0:2] = Te[3:5, 3:5] =  [[n[0], n[1]], [-n[1], n[0]]]
        # Column stiffness
        EA, EI = self.material.E * self.A, self.material.E * self.Izz
        K1 = EA / h * array([[ 1., 0., 0.,-1., 0., 0.],
                             [ 0., 0., 0., 0., 0., 0.],
                             [ 0., 0., 0., 0., 0., 0.],
                             [-1., 0., 0., 1., 0., 0.],
                             [ 0., 0., 0., 0., 0., 0.],
                             [ 0., 0., 0., 0., 0., 0.]])
        # Beam stiffness
        K2 = 2. * EI / h**3 * array([[0.,  0.,    0.,     0.,  0.,   0     ],
                                     [0.,  6.,    3.*h,   0., -6.,   3.*h  ],
                                     [0.,  3.*h,  2.*h*h, 0., -3.*h, h*h   ],
                                     [0.,  0.,    0.,     0.,  0.,   0.    ],
                                     [0., -6.,   -3.*h,   0.,  6.,  -3.*h  ],
                                     [0.,  3.*h,  h*h,    0., -3.*h, 2.*h*h]])
        return dot(dot(Te.T, K1+K2), Te)

# --------------------------------------------------------------------------- #
# ------------------------ PLANE 3 NODE ELEMENT ----------------------------- #
# --------------------------------------------------------------------------- #
class Tria3(Element):
    numdim, numnod = 2, 3
    edges = array([[0,1], [1,2], [2,0]])
    gaussp = array([[1., 1.], [4., 1.], [1., 4.]]) / 6.
    gaussw = ones(3) / 3.
    def isop_map(self, xi):
        xc = self.xc
        x = xc[0,0]*xi[0] + xc[1,0]*xi[1] + xc[2,0]*(1-xi[0]-xi[1])
        y = xc[0,1]*xi[0] + xc[1,1]*xi[1] + xc[2,1]*(1-xi[0]-xi[1])
        return x, y

    def jacobian(self):
        # Element Jacobian
        ((X1, Y1), (X2, Y2), (X3, Y3)) = self.xc
        return (X1*Y2 - X2*Y1 - X1*Y3 + X3*Y1 + X2*Y3 - X3*Y2)

    def shape(self, xp):
        Je = self.jacobian()
        ((X1, Y1), (X2, Y2), (X3, Y3)) = self.xc
        x, y = xp
        # Element shape functions
        Ne = array([X2*Y3 - X3*Y2 + (Y2 - Y3) * x + (X3 - X2) * y,
                    X3*Y1 - X1*Y3 + (Y3 - Y1) * x + (X1 - X3) * y,
                    X1*Y2 - X2*Y1 + (Y1 - Y2) * x + (X2 - X1) * y]) / Je
        return Ne

    def shapegrad(self):
        Je = self.jacobian()
        ((X1, Y1), (X2, Y2), (X3, Y3)) = self.xc
        # Element shape function gradients
        dN = array([[Y2 - Y3, Y3 - Y1, Y1 - Y2],
                    [X3 - X2, X1 - X3, X2 - X1]]) / Je
        return dN

    def edge_shape(self, edge, xp):
        # ordering of nodes
        xb = self.xc[self.edges[edge]]
        he = sqrt((xb[1,0]-xb[0,0])**2 + (xb[1,1]-xb[0,1])**2)
        o = array({0:[0,1,2],1:[2,0,1],2:[1,2,0]}[edge])
        s = he * (xp + 1) / 2.0
        return array([(he - s) / he, s / he, 0.])[o]

class Tria6(Element):
    numdim, numnod = 2, 6
    edges = array([[0,1,3], [1,2,4], [2,0,5]])

class Quad4(Element):
    numdim, numnod = 2, 4
    edges = array([[0,1], [1,2], [2,3], [3,0]])

class Quad8(Element):
    numdim, numnod = 2, 8
    edges = array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]])

# --------------------------------------------------------------------------- #
# ------------------------ HEAT TRANSFER ELEMENT ---------------------------- #
# --------------------------------------------------------------------------- #
class DiffussiveHeatTransfer2D3(Tria3):
    ndof = 1
    signature = 0000001  # 3 NODE 2D HEAT TRANSFER
    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.nodes = asarray(elenod, dtype=int)
        self.xc = asarray(elecoord, dtype=float)
        self.material = elemat
        if elefab:
            raise ValueError('Element does not take element fabrication properties')

    def stiffness(self, *args):
        sfilm = args[0]
        # Material contribution
        Ke = self.stiffness1()
        # Convection contribution
        for iedge in range(len(self.edges)):
            Too, h = sfilm[iedge]
            Ke += self.stiffness2(iedge, h)
        return Ke

    def stiffness1(self):
        # Material stiffness - "resistance" to conduction
        Je = self.jacobian()
        B = self.shapegrad()
        N = [self.shape(self.isop_map(xi)) for xi in self.gaussp]
        w, k = self.gaussw, self.material.isotropic_thermal_conductivity(2)
        return Je/2.*sum([w[i]*dot(dot(B.T,k),B) for (i, Ni) in enumerate(N)], 0)

    def stiffness2(self, edge, h):
        # convection stiffness
        xi, w = array([-1., 1.])/sqrt(3.), ones(2)
        # Determine edge length
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

    def force(self, f, sflux, sfilm):
        Fe = self.heat_source(f)
        for iedge in range(len(self.edges)):
            qn = sflux[iedge]
            if qn:
                # evaluate the boundary flux contribution
                Fe += self.conduction_flux_array(iedge, qn)
            Too, h = sfilm[iedge]
            if h:
                # evaluate the convection contribution
                Fe += self.convection_flux_array(iedge, Too, h)
        return Fe

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

# --------------------------------------------------------------------------- #
# --------------------------- ELEMENT TESTS --------------------------------- #
# --------------------------------------------------------------------------- #
def test_1():
    """Link2 stiffness test"""
    class mat: E = 1
    El = Link1D2(1, [0, 1], [0, 1], mat, A=1)
    K1D = El.stiffness()
    assert allclose([[1,-1],[-1,1]], K1D)

def test_2():
    """Link2 stiffness test"""
    class mat: E=1000
    El = Link2D2(1, [0, 1], [[0,0], [30,40]], mat, A=5)
    K2D = El.stiffness()
    assert allclose([[ 36.,  48., -36., -48.],
                     [ 48.,  64., -48., -64.],
                     [-36., -48.,  36.,  48.],
                     [-48., -64.,  48.,  64.]], K2D)

def test_3():
    """Link2 stiffness test"""
    class mat: E = 343
    El = Link3D2(1, [0, 1], [[0,0,0],[2,3,6]], mat, A=10)
    K3D = El.stiffness()
    assert allclose([[  40.,   60.,  120.,  -40.,  -60., -120.],
                     [  60.,   90.,  180.,  -60.,  -90., -180.],
                     [ 120.,  180.,  360., -120., -180., -360.],
                     [ -40.,  -60., -120.,   40.,   60.,  120.],
                     [ -60.,  -90., -180.,   60.,   90.,  180.],
                     [-120., -180., -360.,  120.,  180.,  360.]], K3D)

def test_4():
    """Beam-Column stiffness test"""
    coord = array([[0, 0], [3, 4]], dtype=float)
    class mat: E = 100
    A, Izz = 125, 250
    El = BeamColumn2D(1, [0, 1], coord, mat, A=A, Izz=Izz)
    Ke = El.stiffness()

if __name__ == '__main__':
    test_1()
    test_2()
    test_3()
    test_4()
