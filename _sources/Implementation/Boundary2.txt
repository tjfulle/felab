Application of Boundary Conditions
==================================

Overview
--------

After assembly, the global system of equations :math:`[K]\{u\}=\{F\}`
that describe the finite element problem is singular. Meaning that, the
equations for the unknown degrees of freedom are non-invertible. This
situation arises when there are no solutions or infinite solutions to a
problem. Boundary conditions must be applied to make the system
non-singular.

Boundary conditions are applied by

-  modifying the global stiffness and force arrays on Dirichlet
   boundaries, and

-  modifying the global force on Neummann boundaries.

The method shown adopts the same assumptions as the assembler in :ref:`Assembly1`.

Example
-------

Consider the problem of determining the displacements along the tapered
bar shown

.. figure:: DiscretizedTaperedBar.png
   :align: center

   Tapered bar discretized by four spring elements.

The assembled global equations are given by

.. math::

     \begin{bmatrix}
       k_{1} & -k_{1} & 0 & 0 & \\
       -k_{1} & k_{1} + k_{2} & -k_{2} & 0 & 0 \\
       0 & -k_{2} & k_{2} + k_{3} & -k_{3} & 0 \\
       0 & 0 & -k_{3} & k_{3} + k_{4} & -k_{4} \\
       0 & 0 & 0 & -k_{4} & k_{4}
     \end{bmatrix}
     \begin{Bmatrix} u_1 \\ u_{2} \\ u_{3} \\ u_{4} \\ u_{5} \end{Bmatrix} =
     \begin{Bmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ \end{Bmatrix}

where :math:`k_e=A_eE_e/h_e` is the stiffness coefficient of the
e\ :math:`^{\rm th}` element. The boundary conditons are

.. math::

   u_{x=0}=\overline{u}_1, \quad \left[AE\frac{du}{dx}\right]_{x=L}=F

Substituting in to gives

.. math::

     \begin{bmatrix}
       k_{1} & -k_{1} & 0 & 0 & \\
       -k_{1} & k_{1} + k_{2} & -k_{2} & 0 & 0 \\
       0 & -k_{2} & k_{2} + k_{3} & -k_{3} & 0 \\
       0 & 0 & -k_{3} & k_{3} + k_{4} & -k_{4} \\
       0 & 0 & 0 & -k_{4} & k_{4}
     \end{bmatrix}
     \begin{Bmatrix}
       \overline{u}_1 \\ u_{2} \\ u_{3} \\ u_{4} \\ u_{5}
     \end{Bmatrix} =
     \begin{Bmatrix}
       R \\ 0 \\ 0 \\ 0 \\ F \\
     \end{Bmatrix}

where :math:`R` is an unkown reaction at node 1. The displacements
:math:`u_{2}, \ u_{3},` and :math:`u_{4}` do not explicitly depend on
:math:`R`, allowing to be written as

.. math::

     \begin{bmatrix}
       1 & 0 & 0 & 0 & 0 \\
       -k_{1} & k_{1} + k_{2} & -k_{2} & 0 & 0 \\
       0 & -k_{2} & k_{2} + k_{3} & -k_{3} & 0 \\
       0 & 0 & -k_{3} & k_{3} + k_{4} & -k_{4} \\
       0 & 0 & 0 & -k_{4} & k_{4}
     \end{bmatrix}
     \begin{Bmatrix}
       u_1 \\ u_{2} \\ u_{3} \\ u_{4} \\ u_{5}
     \end{Bmatrix} =
     \begin{Bmatrix}
       \overline{u}_1 \\ 0 \\ 0 \\ 0 \\ F \\
     \end{Bmatrix}

The global system can now be solved for the unknown displacements and the
unknown reaction can then be determined as :math:`R=-k_1\overline{u}_1 +
\left(k_1+k_2\right)u_2`. However, solving has the disadvantange that it is no
longer a symmetric system due to the way the boundary conditions were applied.
Symmetry of the global equation can be easily recovered by adding
:math:`\begin{bmatrix}0 & k_1\overline{u}_1 & 0 & 0 & 0\end{bmatrix}^T` to
both sides of , giving

.. math::

     \begin{bmatrix}
       1 & 0 & 0 & 0 & 0 \\
       0 & k_{1} + k_{2} & -k_{2} & 0 & 0 \\
       0 & -k_{2} & k_{2} + k_{3} & -k_{3} & 0 \\
       0 & 0 & -k_{3} & k_{3} + k_{4} & -k_{4} \\
       0 & 0 & 0 & -k_{4} & k_{4}
     \end{bmatrix}
     \begin{Bmatrix}
       u_1 \\ u_{2} \\ u_{3} \\ u_{4} \\ u_{5}
     \end{Bmatrix} =
     \begin{Bmatrix}
       \overline{u}_1 \\ k_{1}\overline{u}_1 \\ 0 \\ 0 \\ F \\
     \end{Bmatrix}

In general, for known DOF :math:`\overline{u}_i`, the global stiffness
and force arrays are modified according to

.. math::

     K_{ij} = K_{ji} = \delta_{ij}, \quad
     F_j = \begin{cases}
       \overline{u}_i & i = j \\ F_j-K_{ji}\overline{u}_i & i \ne j
     \end{cases}

where :math:`\delta_{ij}` is the Kronecker delta.

Computational Implementation
----------------------------

In ``pyfem2``, application of the boundary conditions is performed in the method
:ref:`FiniteElementModel.apply_bc_udof <apply_bc>`. (``udof`` stands for
"uniform degree of freedom").

A representative python code that applies boundary conditions is the function
``ApplyBoundaryConditions``, shown below. ``ApplyBoundaryConditions`` is invoked as

.. code:: python

   Kbc, Fbc = ApplyBoundaryConditions(K, F, doftags, dofvals)

The arguments ``ApplyBoundaryConditions`` have been described in :ref:`Assembly1`.

The output is:

+---------+------------------------------------------------+
| ``Kbc`` | Boundary condition modified global stiffness   |
+---------+------------------------------------------------+
| ``Fbc`` | Boundary condition modified global force       |
+---------+------------------------------------------------+

.. code:: python

    def ApplyBoundaryConditions(K, F, doftags, dofvals):

        N, m = doftags.shape
        Kbc, Fbc = K.copy(), F.copy()

        # Dirichlet boundary conditions
        for i in range(N):
            for j in range(m):
                if doftags[i,j] == 0:
                    I = i * m + j
                    Fbc -= [K[k,I] * dofvals[i,j] for k in range(N*m)]
                    Kbc[I, :] = Kbc[:, I] = 0
                    Kbc[I, I] = 1

        # Further modify RHS for Dirichlet boundary
        # This must be done after the loop above.
        for i in range(N):
            for j in range(m):
                if doftags[i,j] == 0:
                    Fbc[i*m+j] = dofvals[i,j]

        return Kbc, Fbc
