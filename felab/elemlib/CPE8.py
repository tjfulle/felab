from .CMDN import CMDN
from .gauss_rule_info import quad_gauss_rule_info


class CPE8(CMDN):
    """8 node plane-strain stress-displacement element

    Notes
    -----
    Node and element face numbering

               [2]
            3---6---2
            |       |
       [3]  7       5 [1]
            |       |
            0---4---1
               [0]

    """
    nodes = 8
    dimensions = 2

    #: The element signature is the active degrees of freedom at each node. Each
    # node signature takes the following form:
    #
    #   (X, Y, Z, TX, TY, TZ, T)
    #
    # ie, the X is the X dof, TX the rotation about the X axis, T the
    # temperature etc. Indicate that the node activates a DOF by assigning a
    # value of 1 to the appropriate slot in the node's signature and 0
    # otherwise.
    #
    #  For a 8 node isoparametric stress-displacement element, the x and y
    #  displacements are the only dofs that are active.
    signature = None

    #: ndir is the number of direct components in the stress tensor. The direct
    # components are the XX, YY, and ZZ components of the stress.
    ndir = None

    #: nshr is the number of shear components in the stress tensor.
    nshr = None

    #: num_gauss is the number of integration points.
    num_gauss = 9

    #: cp defines the center of the element *in the natural coordinates*
    cp = None

    #: xp defines the corner nodes *in the natural coordinates*
    xp = None

    #: edges defines the nodes that define the edges of the element
    edges = None

    def __init__(self, label, nodes, coord, material, t=1.0):
        """8 node plane-strain stress-displacement element

        Parameters
        ----------
        label : int
            The element's label (ID)
        nodes : ndarray
            The node labels making up this element
        coord : ndarray
            The coordinates of the nodes making up this element
        material : Material
            The material model
        t : float
            The plane strain thickness

        """
        super(CPE8, self).__init__(label, nodes, coord, material)
        self.t = t

        # This element must be implemented
        raise NotImplementedError

    @staticmethod
    def gauss_rule_info(point=None):
        return quad_gauss_rule_info(rule=3, point=point)

    @property
    def area(self):
        """Returns the area of the undeformed element"""
        x, y = self.xc[:, [0, 1]].T
        A2 = (x[0] * y[1] - x[1] * y[0]) + (x[1] * y[2] - x[2] * y[1])
        A2 += (x[2] * y[3] - x[3] * y[2]) + (x[3] * y[0] - x[0] * y[3])
        return A2 / 2.0

    @property
    def volume(self):
        """Returns the volume of the undeformed element"""
        return self.t * self.area

    def shape(self, qcoord, edge=None):
        """Evaluate the shape function

        Parameters
        ----------
        qcoord : ndarray
            The natural coordinate
        edge : int
            The edge number, if evaluating the shape function on an edge

        Returns
        -------
        N : ndarray
            The shape function evaluated at qcoord

        """
        raise NotImplementedError

    def shapegrad(self, qcoord):
        """Evaluate the derivative of the shape function

        Parameters
        ----------
        qcoord : ndarray
            The natural coordinate

        Returns
        -------
        dN : ndarray
            The shape function derivative evaluated at qcoord

        """
        raise NotImplementedError

    def shapefun_der(self, coord, qcoord):
        """Shape function and derivative of 8 node quadratic element

        Parameters
        ----------
        coord : ndarray
            The coordinate in the physical coordinates
        qcoord : ndarray
            The coordinate in the natural coordinates

        Returns
        -------
        N : ndarray
            The shape function evaluated at Gauss coordinate
        dNdx : ndarray
            The shape function derivative in the physical coordinates
        J : float
            The Jacobian of the transformation

        """
        raise NotImplementedError

    def bmatrix(self, dN, *args):
        """Return the element's 'B' matrix

        Parameters
        ----------
        dN : ndarray
            The shape function derivative, in the physical coordinates

        Returns
        -------
        B : ndarray
            The element B matrix

        """
        raise NotImplementedError

    def stiffness(self, svars, u, du, time, dtime, kstep, kframe, predef):
        """Evaluate the element stiffness

        Parameters
        ----------
        svars : ndarray
            State variables (stress, strain, etc)
        u : ndarray
            The displacement at the beginning of the step
        du : ndarray
            The displacment increment
        time : float
            The time at the beginning of the step
        dtime : float
            The time increment
        kstep, kframe : int
            The step and frame numbers
        predef : ndarray
            Prescribed deformation fields

        Returns
        -------
        A : ndarray
            The element stiffness

        """
        raise NotImplementedError
