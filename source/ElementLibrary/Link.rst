
.. _BarElements:

Elastic bar elements
====================

+---------------------+---------------------------------+
| ``ElasticLink1D2``  | 1D 2 node elastic bar element   |
+---------------------+---------------------------------+
| ``ElasticLink2D2``  | 2D 2 node elastic bar element   |
+---------------------+---------------------------------+
| ``ElasticLink3D2``  | 3D 2 node elastic bar element   |
+---------------------+---------------------------------+

Active degrees of freedom
-------------------------

-  1D: ``X``

-  2D: ``X``, ``Y``

-  3D: ``X``, ``Y``, ``Z``

Element formulation
-------------------

The two-node elastic bar element is shown embedded in a two-dimensional
space in Figure [fig:ch1.link]. In the element’s local coordinates
(:math:`\tilde{\pmb{x}}_i`), forces and displacements act along only the
local :math:`\tilde{x}` direction. When viewed in the global
coordinates, forces and displacements act along all :math:`\pmb{x}_i`
coordinate directions. Thus, in two dimensions, the bar element has
:math:`2` degrees of freedom per node, for a total of :math:`4` degrees
of freedom. Elements of this type are referred to as **plane** bar
elements and a system made up of plane bar elements is a **plane
truss**. In three dimensions, the element has :math:`3` degrees of
freedom per node, for a total of 6 degrees of freedom and is referred to
as **space** bar elements. A system made up of space bar elements is a
**space truss**.

.. figure:: TrussCoords.png
   :align: center

   Two-node elastic bar element

In three dimensions, the element node displacement and force vectors in
three dimensions are arranged as

.. math::

   \boldsymbol{u}-{(e)} = \begin{Bmatrix}
       u_{1x} \\ u_{1y} \\ u_{1z} \\ u_{2x} \\ u_{2y} \\ u_{2z}
     \end{Bmatrix}, \quad
     \boldsymbol{f}-{(e)} = \begin{Bmatrix}
       f_{1x} \\ f_{1y} \\ f_{1z} \\ f_{2x} \\ f_{2y} \\ f_{2z}
     \end{Bmatrix}

The displacement and force vectors in the global coordinates are related
to their local counterparts through the following transformations

.. math::

   \tilde{\boldsymbol{u}}-{(e)} = \boldsymbol{T}\boldsymbol{u}-{(e)},
     \quad
     \tilde{\boldsymbol{f}}-{(e)} = \boldsymbol{T}\boldsymbol{f}-{(e)},

where :math:`\boldsymbol{T}` is the :math:`2 \times 6` matrix of
direction cosines of the axial axis of the element, given by

.. math::

   \label{eq:transm}
     \boldsymbol{T} = \begin{bmatrix}
       l_{ij} & m_{ij} & n_{ij} & 0 & 0 & 0 \\
       0 & 0 & 0 & l_{ij} & m_{ij} & n_{ij} \\
     \end{bmatrix}, \quad

The direction cosines :math:`l_{ij}`, :math:`m_{ij}`, and :math:`n_{ij}`
are given by

.. math::

   \begin{aligned}
     \label{eq:dircos}
     l_{ij} &= \cos\left(\tilde{x}, x\right) = \frac{x_j-x_i}{l_e} \\
     m_{ij} &= \cos\left(\tilde{x}, y\right) = \frac{y_j-y_i}{l_e} \\
     n_{ij} &= \cos\left(\tilde{x}, y\right) = \frac{z_j-z_i}{l_e}\end{aligned}

The force-displacement relationship in the local coordinates is

.. math::

   \tilde{\boldsymbol{f}}-{(e)} =
     \tilde{\boldsymbol{k}}-{(e)}\tilde{\boldsymbol{u}}-{(e)}

where :math:`\tilde{\boldsymbol{k}}-{(e)}` is the :math:`2 \times 2`
element stiffness matrix in the local coordinates. Using the coordinate
transformations, the force-displacement relationship in the global
coordinates is

.. math::

   \label{eq:ch1.force-disp-1}
     \tilde{\boldsymbol{f}}-{(e)}
     = \boldsymbol{T}\boldsymbol{f}-{(e)}
     = \tilde{\boldsymbol{k}}-{(e)}\tilde{\boldsymbol{u}}-{(e)}
     = \tilde{\boldsymbol{k}}-{(e)}\boldsymbol{T}\boldsymbol{u}-{(e)}

From , we get

.. math:: \boldsymbol{f}-{(e)} = \boldsymbol{k}-{(e)}\boldsymbol{u}-{(e)}

where the element stiffness matrix in the global coordinates
:math:`\boldsymbol{k}-{(e)}` is

.. math::

   \label{eq:gstiff1}
     \boldsymbol{k}-{(e)} = \boldsymbol{T}-T\tilde{\boldsymbol{k}}-{(e)}\boldsymbol{T}

Expanding in components, the force-displacement relationship can be
expressed as

.. math::

   \begin{Bmatrix}
       \boldsymbol{f}_1 \\ \boldsymbol{f}_2
     \end{Bmatrix} = \frac{AE}{l_e}
     \begin{bmatrix*}[r]
       \boldsymbol{n}\otimes\boldsymbol{n} &
       -\boldsymbol{n}\otimes\boldsymbol{n} \\
       -\boldsymbol{n}\otimes\boldsymbol{n} &
       \boldsymbol{n}\otimes\boldsymbol{n} \\
     \end{bmatrix*}
     \begin{Bmatrix}
       \boldsymbol{u}_1 \\ \boldsymbol{u}_2
     \end{Bmatrix}

where the vectors :math:`\boldsymbol{f}_i` and :math:`\boldsymbol{u}_i`
are the force and displacement vectors for the :math:`i-{\rm th}` node
in the global coordinates, :math:`A` and :math:`E` are the constant
element area and Young’s modulus, respectively, :math:`l_e` is the
element length, and :math:`\boldsymbol{n}` is the element unit normal,
given by

.. math::

   \begin{aligned}
     \boldsymbol{v} &= \boldsymbol{x}_2-\boldsymbol{x}_1 \\
     \boldsymbol{n} &= \frac{\boldsymbol{v}}{\lVert\boldsymbol{v}\rVert}\end{aligned}

The product :math:`\boldsymbol{a}\otimes\boldsymbol{b}` represents the
vector outer product whose result is a matrix. In direct, indicial, and
matrix notation, the matrix :math:`\boldsymbol{A}` resulting from the
outer product of vectors :math:`\boldsymbol{a}` and
:math:`\boldsymbol{b}` is

.. math::

   \boldsymbol{A} = \boldsymbol{a}\otimes\boldsymbol{b} = a_i b_j =
     \begin{bmatrix}
       a_1b_1 & a_1b_2 & \ldots & a_1b_n \\
       a_2b_1 & a_2b_2 & \ldots & a_2b_n \\
       \vdots & \vdots & \ddots & \vdots \\
       a_nb_1 & a_nb_2 & \ldots & a_nb_n
     \end{bmatrix}

Element stiffness function
--------------------------

The computation of stiffness matrix for the two-node elastic bar is
performed in the function , given in Listing [lst:stiff1]. is invoked as

.. code:: python

    ke = Link2Stiffness(xc, E, A)

Arguments to the are

+----------+---------------------+
| ``xc``   | Nodal coordinates   |
+----------+---------------------+
| ``E``    | Element modulus.    |
+----------+---------------------+
| ``A``    | Element area.       |
+----------+---------------------+

The output from is

+----------+-------------------------------------------------------------------------------------------------------------------------------------+
| ``ke``   | Element stiffness stored as a (numdim\*2, numdim\*2) symmetric matrix, where numdim is the number of degrees of freedom per node.   |
+----------+-------------------------------------------------------------------------------------------------------------------------------------+

.. code:: python

    def Link2Stiffness(xc, E, A):
        # Element dimensionality
        xc = asarray(xc)
        if xc.ndim == 1:
            numdim = 1
        else:
            numdim = xc.shape[1]
        # Compute element normal
        v = xc[1] - xc[0]
        h = sqrt(dot(v, v))
        n = v / h
        if xc.ndim == 1:
            nn = 1.
        else:
            nn = outer(n, n)
        # Assemble element stiffness
        k = zeros((2*numdim, 2*numdim))
        i, j = numdim, 2*numdim
        k = zeros((2*numdim, 2*numdim))
        k[0:i, 0:i] = k[i:j, i:j] =  nn # upper left and lower right 2x2
        k[0:i, i:j] = k[i:j, 0:i] = -nn # lower left and upper right 2x2
        return A * E / h * k

Verification of Program Units
-----------------------------

The validity of each program unit is tested below. In each case, the
Python statement is used with the function to test the output. The
statement performs a Null operation if the expression that follows
evaluates to , otherwise it raises an error. The function evaluates to
if all elements in two test arrays are close to within a toloerance,
otherwise it evaluates to .

Verification of the n-Dimensional Elastic Bar Stiffness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    K1D = Link2Stiffness([0, 1], 1, 1)
    assert allclose([[1,-1],[-1,1]], K1D)

.. code:: python

    K2D = Link2Stiffness([[0,0], [30,40]], 5, 1000)
    assert allclose([[ 36.,  48., -36., -48.],
                     [ 48.,  64., -48., -64.],
                     [-36., -48.,  36.,  48.],
                     [-48., -64.,  48.,  64.]], K2D)

.. code:: python

    K3D = Link2Stiffness([[0,0,0],[2,3,6]], 10, 343)
    assert allclose([[  40.,   60.,  120.,  -40.,  -60., -120.],
                     [  60.,   90.,  180.,  -60.,  -90., -180.],
                     [ 120.,  180.,  360., -120., -180., -360.],
                     [ -40.,  -60., -120.,   40.,   60.,  120.],
                     [ -60.,  -90., -180.,   60.,   90.,  180.],
                     [-120., -180., -360.,  120.,  180.,  360.]], K3D)
