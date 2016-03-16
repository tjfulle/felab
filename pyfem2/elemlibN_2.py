from numpy import *
from .utilities import *
from .elemlib1 import Element

__all__ = ['ElasticLink1D2', 'ElasticLink2D2', 'ElasticLink3D2', 'BeamColumn2D']

# --------------------------------------------------------------------------- #
# ------------------------------ TRUSS ELEMENT ------------------------------ #
# --------------------------------------------------------------------------- #
class ElasticLinknD2(Element):
    """Base class for 2 node elastic link elements

    Parameters
    ----------
    label : int
        Element label
    elenod : list of int
        Internal node IDs of nodes making up this element
    elecoord : ndarray
        Coordinates of nodes
    elemat : object
        A pyfem2.mat.Material instance
    elefab : dict
        Requires area 'A'

    """
    ndir, nshr = 1, 0
    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.inodes = asarray(elenod, dtype=int)
        self.xc = asarray(elecoord, dtype=float)
        self.material = elemat
        self.A = elefab.get('A')
        if self.A is None:
            raise UserInputError("Expected exactly area 'A' as the only element "
                                 "fabrication property")

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type, load_fac):
        """Computes the response of a n-dimensional elastic link

        Parameters
        ----------

        Returns
        -------
        k : array_like
            (2*nd, 2*nd) elastic stiffness
        f : array_like
            (2*nd, 0) internal force

        """
        # INTERNAL FORCE
        ndof = count_digits(self.signature[0])

        compute_stiff = cflag in (STIFF_AND_FORCE, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_FORCE, FORCE_ONLY)

        if compute_force:
            Fe = zeros(2*ndof)

        if cflag == FORCE_ONLY:
            return Fe

        # ELEMENT NORMAL
        v = self.xc[1] - self.xc[0]
        h = sqrt(dot(v, v))
        n = v / h
        if self.dimensions == 1:
            nn = 1.
        else:
            nn = outer(n, n)

        # ASSEMBLE ELEMENT STIFFNESS
        i, j = ndof, 2*ndof
        Ke = zeros((2*ndof, 2*ndof))
        Ke[0:i, 0:i] = Ke[i:j, i:j] =  nn # UPPER LEFT AND LOWER RIGHT 2X2
        Ke[0:i, i:j] = Ke[i:j, 0:i] = -nn # LOWER LEFT AND UPPER RIGHT 2X2
        Ke *= self.A * self.material.E / h

        if cflag == STIFF_AND_FORCE:
            return Ke, Fe

        elif cflag == STIFF_ONLY:
            return Ke

    def internal_force(self, uc):
        """
        .. _link_int_force:

        Computes the element axial internal force

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

class ElasticLink1D2(ElasticLinknD2):
    nodes = 2
    variables = ('P', 'S')
    dimensions = 1
    signature = [(1,0,0,0,0,0,0),  # 2 NODE 1D LINE LINK
                 (1,0,0,0,0,0,0)]

class ElasticLink2D2(ElasticLinknD2):
    nodes = 2
    variables = ('P', 'S')
    dimensions = 2
    signature = [(1,1,0,0,0,0,0),  # 2 NODE 2D LINE LINK
                 (1,1,0,0,0,0,0)]

class ElasticLink3D2(ElasticLinknD2):
    nodes = 2
    variables = ('P', 'S')
    dimensions = 3
    signature = [(1,1,1,0,0,0,0),
                 (1,1,1,0,0,0,0)]  # 2 NODE 3D LINE LINK

# --------------------------------------------------------------------------- #
# ---------------------- EULER BERNOULI BEAM ELEMENT ------------------------ #
# --------------------------------------------------------------------------- #
class BeamColumn2D(Element):
    """Base class for 2 node Euler-Bernouli beam elements

    Parameters
    ----------
    label : int
        Element label
    elenod : list of int
        Internal node IDs of nodes making up this element
    elecoord : ndarray
        Coordinates of nodes
    elemat : object
        A pyfem2.mat.Material instance
    elefab : dict
        Requires area 'A' and 'Izz'

    """
    ndir, nshr = 1, 0
    nodes = 2
    variables = ('P', 'S')
    dimensions = 2
    signature = [(1,1,0,0,0,1,0),
                 (1,1,0,0,0,1,0)]
    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.inodes = asarray(elenod, dtype=int)
        self.xc = asarray(elecoord, dtype=float)
        self.material = elemat
        self.A, self.Izz = elefab.get('A'), elefab.get('Izz')
        if self.A is None or self.Izz is None:
            raise ValueError('Incorrect element fabrication properties')

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type, load_fac):

        # INTERNAL FORCE
        Fe = zeros(6)

        if cflag == FORCE_ONLY:
            return Fe

        # COMPUTE ELEMENT NORMAL
        v = self.xc[1] - self.xc[0]
        h = sqrt(dot(v, v))
        n = v / h

        # TRANSFORMATION MATRIX
        Te = eye(6)
        Te[0:2, 0:2] = Te[3:5, 3:5] =  [[n[0], n[1]], [-n[1], n[0]]]

        # COLUMN STIFFNESS
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

        Ke = dot(dot(Te.T, K1+K2), Te)

        if cflag == STIFF_AND_FORCE:
            return Ke, Fe

        elif cflag == STIFF_ONLY:
            return Ke
