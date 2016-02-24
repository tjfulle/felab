from numpy import *
from numpy.linalg import solve

def UniformBar(xa, xb, A, E, p, q, numele=10):
    """Solve for the displacement along a uniform bar, fixed on its left end
    and with a constant distributed load acting over its length and point
    force applied to its right, using the finite element method.

    Parameters
    ----------
    xa, xb : float
        Left and right coordinates of bar
    A, E : float
        Area and Young's modulus of bar
    p : float
        Point force acting on right end of bar
    q : float
        Distributed load acting over bar
    numele : int, optional {10}
        Number of elements

    Returns
    -------
    u : array_like
        Displacements of nodes along bar
    e : array_like
        Element strains
    f : array_like
        Element forces

    Examples
    --------

    Define a uniform bar on the interval :math:`x\in (0, 1)` with cross
    sectional area :math:`A=1`, Young's modulus :math:`E=1`, point force
    :math:`p=1`, and distributed load :math:`q=0`:

    >>> xa, xb = 0., 1.
    >>> A, E = 1., 1.
    >>> p, q = 1., 0.
    >>> u, e, f = UniformBar(xa, xb, A, E, p, q, numele=10)
    >>> assert abs(u[-1] - 1.) < 1e-12

    """

    # Number of nodes, coordinates, and element connectivity
    numnod = numele + 1
    coord = linspace(xa, xb, numnod)
    elecon = array([[el, el+1] for el in range(numele)])

    # loop over elements in the connectivity and assemble global stiffness, force
    K, F = zeros((numnod, numnod)), zeros(numnod)
    for iel in range(numele):
        # Element length, stiffness, and force
        nodes = elecon[iel]
        xe = coord[nodes]
        he = xe[1] - xe[0]
        ke = A * E / he * array([[1., -1.], [-1., 1.]])
        qe = q * he / 2. * ones(2)  # distributed load contribution

        # add contributions to global matrices
        F[nodes] += qe
        for i in range(2):
            I = nodes[i]
            for j in range(i, 2):
                J = nodes[j]
                K[I,J] += ke[i,j]
                K[J,I] = K[I,J]

    # Apply boundary conditions
    K[0, :] = 0.
    K[0, 0] = 1.
    F[0] = 0
    F[-1] = p

    # Solve for u
    u = solve(K, F)

    # determine element forces and strains
    elefor, e = zeros(numele), zeros(numele)
    for i in range(numele):
        du = u[i+1] - u[i]
        he = coord[i+1] - coord[i]
        elefor[i] = A * E / he * du
        e[i] = du / he

    return u, e, elefor
