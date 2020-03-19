from .CMDN import CMDN
from .gauss_rule_info import tri_gauss_rule_info


class CPE3(CMDN):
    """3-node isoparametric plane strain stress-displacement element

    Notes
    -----
    Node and element face numbering


            1
            | .
       [1]  |   .  [0]
            |     .
            2------0
              [2]

    """
    nodes = 3
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
    signature = None

    #: ndir is the number of direct components in the stress tensor. The direct
    # components are the XX, YY, and ZZ components of the stress. For plane
    # strain, there are 3 direct components.
    ndir = None

    #: nshr is the number of shear components in the stress tensor. Only the xy
    # component is nonzero in plane strain, so nshr is 1
    nshr = None

    #: num_gauss is the number of integration points. The CPE4 is fully
    # integrated and, thus, has 4 integration points
    num_gauss = None

    #: cp defines the center of the element *in the natural coordinates*
    cp = None

    #: xp defines the corner nodes *in the natural coordinates*
    xp = None

    #: edges defines the nodes that define the edges of the element
    edges = None

    def __init__(self, label, nodes, coord, material, t=1.0):
        """Initialize the element

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
        super(CPE3, self).__init__(label, nodes, coord, material)
        self.t = t

    @staticmethod
    def gauss_rule_info(point=None):
        """Evaluates the Gauss quadrature rule for this point

        Parameters
        ----------
        point : int
            The Gauss point number

        Returns
        -------
        xi : ndarray
            The natural coordinate associated with this point
        w : float
            The integration weight associated with this point

        """
        return tri_gauss_rule_info(3, point)

    @property
    def area(self):
        """Returns the area of the undeformed element"""
        x, y = self.xc[:, [0, 1]].T
        a = 0.5 * (x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]))
        return a

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
        """Shape function and derivative of 4 node bilinear element

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
