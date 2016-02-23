.. _FirstProgram:

A First Finite Element Program
==============================

Overview
--------

A finite element program for computing the displacements in a
uniform elastic bar fixed at its origin and subject to a uniform body
load and constant point force at its end is developed.

Uniform Bar Program
-------------------

Problem Definition
~~~~~~~~~~~~~~~~~~

Determine the displacements in uniform elastic bar fixed at its origin and
subject to a uniform body load and constant point force at its end, as shown,

.. figure:: FirstTruss.png
   :align: center

   Uniform truss subject to constant body load and point force. Shown is
   the physical model and the 1D finite element idealization. Node
   numbers are shown adjacent to nodes and range from 0 to the number of
   elements :math:`N`. Node numbers are shown circled under their
   element and range from 1 to :math:`N`.

The finite element equations for the elastic bar are:

.. math::

     k_{ij}^eu_j^e = f_i^e

where :math:`k_{ij}^e` is the constant element stiffness, given by

.. math::

     \left[k^e\right] = \frac{A_eE_e}{h_e}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}

and the nodal force :math:`f_i^e` is

.. math::

     \left\{f^e\right\} = \frac{q_0h}{2}\begin{Bmatrix}1\\1\end{Bmatrix}

The assembled global equations are

.. math::

     \begin{bmatrix}
       k_{11}^1 & k_{12}^1 & \\
       k_{21}^1 & k_{22}^1+k_{11}^2 & k_{12}^2 \\
       && \ddots \\
       &&& k_{22}^{i-1}+k_{11}^i & k_{12}^i \\
       &&& k_{21}^{i} & k_{22}^i+k_{11}^{i+1} & k_{12}^{i+1} \\
       &&&& \ddots \\
       &&&&& k_{21}^N & k_{22}^N
     \end{bmatrix}
     \begin{Bmatrix}
       u_0 \\ u_1 \\ \vdots \\ u_i \\ u_{i+1} \\ \vdots \\ u_N
     \end{Bmatrix} =
     \begin{Bmatrix}
       f_0 \\ f_1 \\ \vdots \\ f_i \\ f_{i+1} \\ \vdots \\ f_N
     \end{Bmatrix}

On applying the boundary conditions, the global equations in become

.. math::

     \begin{bmatrix}
       1 & 0 & \\
       k_{21}^1 & k_{22}^1+k_{11}^2 & k_{12}^2 \\
       && \ddots \\
       &&& k_{22}^{i-1}+k_{11}^i & k_{12}^i \\
       &&& k_{21}^{i} & k_{22}^i+k_{11}^{i+1} & k_{12}^{i+1} \\
       &&&& \ddots \\
       &&&&& k_{21}^N & k_{22}^N
     \end{bmatrix}
     \begin{Bmatrix}
       u_0 \\ u_1 \\ \vdots \\ u_i \\ u_{i+1} \\ \vdots \\ u_N
     \end{Bmatrix} =
     \begin{Bmatrix}
       0 \\ f_1 \\ \vdots \\ f_i \\ f_{i+1} \\ \vdots \\ f_N+f
     \end{Bmatrix}

The global system of equations is formed and solved in the program in the
program ``UniformBar``. While simplistic, the program demonstrates the
computational implementation of the preprocessing, processing, and
postprocessing steps outlined above.

The program is invoked as:

.. code:: python

    u, e, elefor = UniformBar(xa, xb, A, E, f, q0, numele)

The arguments to ``UniformBar`` are:

+--------------+---------------------------------------------------------------+
| ``xa``       | The coordinate location of the bar’s origin.                  |
+--------------+---------------------------------------------------------------+
| ``xb``       | The coordinate location of the bar’s end.                     |
+--------------+---------------------------------------------------------------+
| ``A``        | The cross sectional area of the bar.                          |
+--------------+---------------------------------------------------------------+
| ``E``        | The Young’s modulus of the bar.                               |
+--------------+---------------------------------------------------------------+
| ``f``        | The point force applied to the bar’s end.                     |
+--------------+---------------------------------------------------------------+
| ``q0``       | The constant body load applied along the length of the bar.   |
+--------------+---------------------------------------------------------------+
| ``numele``   | Number of elements (optional, default: 10).                   |
+--------------+---------------------------------------------------------------+

The outputs from ``UniformBar`` are:

+--------------+---------------------------+
| ``u``        | The nodal displacements   |
+--------------+---------------------------+
| ``e``        | The element strains       |
+--------------+---------------------------+
| ``elefor``   | The element forces        |
+--------------+---------------------------+

.. code:: python

    from numpy import *
    from numpy.linalg import solve
    def UniformBar(xa, xb, A, E, f, q0, numele=10):
        numnod = numele + 1
        coord = linspace(xa, xb, numnod)
        elecon = array([[el, el+1] for el in range(numele)])
        K, F = zeros((numnod, numnod)), zeros(numnod)

        # loop over elements and assemble global stiffness, force
        for iel in range(numele):
            # Element length, stiffness, and force
            nodes = elecon[iel]
            xe = coord[nodes]
            he = xe[1] - xe[0]
            ke = A * E / he * array([[1., -1.], [-1., 1.]])
            qe = q0 * he / 2. * ones(2)

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
        F[-1] = f

        # Solve for u
        u = solve(K, F)

        # determine element forces and strains
        elefor, e = zeros(numele), zeros(numele)
        for i in range(numele):
            du = u[i+1] - u[i]
            dx = coord[i+1] - coord[i]
            elefor[i] = k * du
            e[i] = du / dx

        return u, e, elefor

How Does it Work?
~~~~~~~~~~~~~~~~~

This section describes each part of ``UniformBar``.

.. code:: python

    from numpy import *
    from numpy.linalg import solve

These statements make the objects and the linear solver accessible to
the program.

.. code:: python

    def UniformBar(xa, xb, A, E, f, q0, numele=10):

This statement defines the function to be a function of and (with the
default being 10).

.. code:: python

    numnod = numele + 1
    coord = linspace(xa, xb, numnod)
    elecon = array([[el, el+1] for el in range(numele)])

These statements declare the number of nodes (one more than the number
of elements), initialize the array of nodal coordinates, and initialize
the element connectivity.

.. code:: python

    K = zeros((numnod, numnod))
    F = zeros(numnod)

These statements initialize the global force and stiffness.

In the following loop, each element’s stiffness and force arrays are
formed and merged in to the global stiffness and force. The merging maps
the appropriate local degrees of freedom to their global counterparts.
Since nodes and elements are ordered continuously, the global DOF is the
local DOF.

.. code:: python

    for iel in range(numele):
        nodes = elecon[iel]
        xe = coord[nodes]
        he = xe[1] - xe[0]
        ke = A * E / he * array([[1., -1.], [-1., 1.]])
        qe = q0 * he / 2. * ones(2)
        F[nodes] += qe
        for i in range(2):
            I = nodes[i]
            for j in range(i, 2):
                J = nodes[j]
                K[I,J] += ke[i,j]
                K[J,I] = K[I,J]

Boundary conditions are applied setting the all but the first components
of the first row of the global stiffness to 0. The first component is
set to 1 and the first component of the global force is set to 0, in
accordance to the displacement boundary condition at the bar’s origin.
Finally, the point force is added to the last component of the global
force.

.. code:: python

    K[0, :] = 0.
    K[0, 0] = 1.
    F[0] = 0.
    F[-1] += f

Nodal displacements are computed.

.. code:: python

    u = solve(K, F)

Finally, the element forces and strains are determined and returned

.. code:: python

    elefor, e = zeros(numele), zeros(numele)
    for i in range(numele):
        du = u[i+1] - u[i]
        dx = coord[i+1] - coord[i]
        elefor[i] = k * du
        e[i] = du / dx
    return u, e, elefor

Conclusion
----------

In the chapters to follow, more sophisticated finite element programs
will be developed. Being more complex not withstanding, they will share
many of the same ideas and patterns of the program outlined above.
