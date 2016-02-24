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
    nfab = 1
    def __init__(self, label, elenod, elecoord, elemat, **elefab):
        self.label = label
        self.nodes = asarray(elenod, dtype=int)
        self.xc = asarray(elecoord, dtype=float)
        self.material = elemat
        self.A = elefab.get('A')
        if self.A is None:
            raise UserInputError("Expected exactly area 'A' as the only element "
                                 "fabrication property")

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
    numdim, ndof, numnod = 1, 1, 2
    signature = (1,0,0,0,0,0,0)  # 2 NODE 1D LINE LINK

class ElasticLink2D2(ElasticLinknD2):
    numdim, ndof, numnod = 2, 2, 2
    signature = (1,1,0,0,0,0,0)  # 2 NODE 2D LINE LINK

class ElasticLink3D2(ElasticLinknD2):
    numdim, ndof, numnod = 3, 3, 2
    signature = (1,1,1,0,0,0,0)  # 2 NODE 3D LINE LINK

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
    nfab = 2
    numdim, ndof, numnod = 2, 3, 2
    signature = (1,1,0,0,0,1,0)

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
