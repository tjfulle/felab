import numpy as np
from .CMDN import CMDN
from .gauss_rule_info import quad_gauss_rule_info


class CPE4(CMDN):
    """4-node isoparametric element stress-displacement plane strain element

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
    nodes = 4
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
    #  For the 4 node isoparametric plane strain element, the x and y
    #  displacements are the only dofs that are active.
    signature = [
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0, 0),
    ]

    #: ndir is the number of direct components in the stress tensor. The direct
    # components are the XX, YY, and ZZ components of the stress. For plane
    # strain, there are 3 direct components.
    ndir = 3

    #: nshr is the number of shear components in the stress tensor. Only the xy
    # component is nonzero in plane strain, so nshr is 1
    nshr = 1

    #: num_gauss is the number of integration points. The CPE4 is fully
    # integrated and, thus, has 4 integration points
    num_gauss = 4

    #: cp defines the center of the element *in the natural coordinates*
    cp = np.array([0, 0], dtype=float)

    #: xp defines the corner nodes *in the natural coordinates*
    xp = np.array(
        [[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0]],
        dtype=float,
    )

    #: edges defines the nodes that define the edges of the element
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])

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
        super(CPE4, self).__init__(label, nodes, coord, material)
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
        return quad_gauss_rule_info(2, point)

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
        if edge is not None:
            # EVALUATE SHAPE FUNCTION ON SPECIFIC EDGE
            xi = qcoord
            qcoord = np.array([[xi, -1.0], [1.0, xi], [xi, 1.0], [-1.0, xi]][edge])
        xi, eta = qcoord
        N = np.array(
            [
                (1.0 - xi) * (1.0 - eta),
                (1.0 + xi) * (1.0 - eta),
                (1.0 + xi) * (1.0 + eta),
                (1.0 - xi) * (1.0 + eta),
            ]
        )
        return N / 4.0

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
        xi, eta = qcoord
        dN = np.array(
            [
                [-1.0 + eta, 1.0 - eta, 1.0 + eta, -1.0 - eta],
                [-1.0 + xi, -1.0 - xi, 1.0 + xi, 1.0 - xi],
            ]
        )
        return dN / 4.0

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

        # Shape function and derivative at Gauss point
        N = self.shape(qcoord)
        dNdxi = self.shapegrad(qcoord)

        # Jacobian to natural coordinates
        dxdxi = np.dot(dNdxi, coord)
        dxidx = np.linalg.inv(dxdxi)
        J = np.linalg.det(dxdxi)

        # Convert shape function derivatives to derivatives wrt global physical
        # coordinates
        dNdx = np.dot(dxidx, dNdxi)

        return N, dNdx, J

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
        B = np.zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B

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
        xc = self.xc

        # Data for indexing state variable array
        ntens = self.ndir + self.nshr
        v = [x[0] for x in self.variables()]
        m = len(v) * ntens
        a1, a2, a3 = [v.index(x) for x in ("E", "DE", "S")]

        # Compute integration point data
        ndof = self.nodes * self.numdofpernod
        A = np.zeros((ndof, ndof))
        for p in range(self.num_gauss):

            # Index to start of state variables
            ij = m * p
            xi, wt = self.gauss_rule_info(p)
            Ne, dNdx, J = self.shapefun_der(xc, xi)

            B = self.bmatrix(dNdx, Ne)

            # STRAIN INCREMENT
            de = np.dot(B, du)

            # predef and increment
            temp = np.dot(Ne, predef[0, 0])
            dtemp = np.dot(Ne, predef[1, 0])

            # Material response
            xv = np.zeros(1)
            e = np.array(svars[0, ij + a1 * ntens : ij + (a1 + 1) * ntens])
            s = np.array(svars[0, ij + a3 * ntens : ij + (a3 + 1) * ntens])
            s, xv, D = self.material.eval(
                s, xv, e, de, time, dtime, temp, dtemp, self.ndir, self.nshr, ntens
            )

            # store the updated variables
            svars[1, ij + a1 * ntens : ij + (a1 + 1) * ntens] += de  # STRAIN
            svars[1, ij + a2 * ntens : ij + (a2 + 1) * ntens] = de  # STRAIN INCREMENT
            svars[1, ij + a3 * ntens : ij + (a3 + 1) * ntens] = s  # STRESS

            # add contribution of function call to integral
            c = J * wt
            A += c * np.dot(np.dot(B.T, D), B)

        return A
