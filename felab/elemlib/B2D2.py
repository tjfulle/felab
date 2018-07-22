from numpy import *

from ..constants import *
from .element_base import element_base

# --------------------------------------------------------------------------- #
# ---------------------- EULER BERNOULI BEAM ELEMENT ------------------------ #
# --------------------------------------------------------------------------- #
class B2D2(element_base):
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
        A felab.mat.Material instance
    elefab : dict
        Requires area 'A' and 'Izz'

    """
    nodes = 2
    dimensions = 2
    signature = [(1,1,0,0,0,1,0),
                 (1,1,0,0,0,1,0)]
    elefab = {'A': None, 'Izz': None}

    @classmethod
    def variables(cls):
        return (('P', SCALAR), ('S', SCALAR))

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type):

        if cflag == LP_OUTPUT:
            return
            raise NotImplementError

        if cflag == RHS_ONLY:
            return zeros(6), zeros(6)

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
        # BEAM STIFFNESS
        K2 = 2. * EI / h**3 * array([[0.,  0.,    0.,     0.,  0.,   0     ],
                                     [0.,  6.,    3.*h,   0., -6.,   3.*h  ],
                                     [0.,  3.*h,  2.*h*h, 0., -3.*h, h*h   ],
                                     [0.,  0.,    0.,     0.,  0.,   0.    ],
                                     [0., -6.,   -3.*h,   0.,  6.,  -3.*h  ],
                                     [0.,  3.*h,  h*h,    0., -3.*h, 2.*h*h]])

        Ke = dot(dot(Te.T, K1+K2), Te)

        if cflag == STIFF_AND_RHS:
            return Ke, zeros(6), zeros(6)

        elif cflag == STIFF_ONLY:
            return Ke
