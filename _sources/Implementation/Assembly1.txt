.. _Assembly1:

Assembly, Part I
================

Overview
--------

Assembly is the process of merging element equations in to their
corresponding locations in the global equations. The assembler described
constructs and merges one element at a time and adopts the following
assumptions:

- internal node numbers are ordered continuously from
  :math:`0\rightarrow n-1`, where :math:`n` is the total number of nodes (though
  external node labels can be assigned arbitrarily, see
  :ref:`InternalNodesAndElements`);

- the number and type of degrees of freedom at each node is the same;

- there are no multifreedom constraints; and

- the global stiffness matrix is stored as a full symmetric matrix.

Assembly by Example
-------------------

The assembly process is clearly described by example.  Consider the following truss:

.. figure:: Truss1.png
   :align: center
   :scale: 150

   Example two element truss.

Each node has two degrees of freedom (translation in the :math:`x` and
:math:`y` directions), giving :math:`3` nodes :math:`\times` 2 DOF/node = 6
DOF for the structure. The array of generalized displacements :math:`\{u\}` is
organized as

.. math::

   \{u\} = \begin{Bmatrix}
       u_x^0 \\ u_y^0 \\ u_x^1 \\ u_y^1 \\ u_x^2 \\ u_y^2
     \end{Bmatrix} =
     \begin{Bmatrix}
       u_0 \\ u_1 \\ u_2 \\ u_3 \\ u_4 \\ u_5
     \end{Bmatrix}

The associated global stiffness :math:`[K]` is a :math:`6\times6` matrix that we
initialize with zeros:

.. math::

   \label{ass.K.1}
     \left[K\right] = \begin{bmatrix}
       0 & 0 & 0 & 0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0
     \end{bmatrix}

Element 0 joins nodes 0 and 1. From the definition of the generalized
displacement array, the local-to-global-DOF mapping for nodes 0 and 1 is:

+--------+-------------+--------------+
| Node   | Local DOF   | Global DOF   |
+========+=============+==============+
| 0      | :math:`x`   | 0            |
+--------+-------------+--------------+
| 0      | :math:`y`   | 1            |
+--------+-------------+--------------+
| 1      | :math:`x`   | 2            |
+--------+-------------+--------------+
| 1      | :math:`y`   | 3            |
+--------+-------------+--------------+

.. note::

   The local to global degree of freedom mapping is stored in the "Element
   freedom table" ``eft``, where ``eft[i]`` is the global degree of freedom
   for local degree of freedom ``i``. For example, the ``eft`` for the
   preceding table is ``eft=[0,1,2,3]``.

The element stiffness matrix :math:`\left[k^0\right]` is a fully populated
matrix given on the left below with the element freedom table entries marking
the columns and rows. The result of merging the element stiffness
:math:`\left[k^0\right]` in to the global stiffness :math:`\left[K\right]` is
shown on the right.

.. math::

   \label{ass.k1.1}
     \stackrel{\begin{matrix}
         \hphantom{k_{11}^0} &
         \hphantom{k_{12}^0} &
         \hphantom{k_{13}^0} &
         \hphantom{k_{14}^0} \\ 0 & 1 & 2 & 3 \end{matrix}}{
     \begin{bmatrix}
       k_{11}^0 & k_{12}^0 & k_{13}^0 & k_{14}^0 \\
       k_{21}^0 & k_{22}^0 & k_{23}^0 & k_{24}^0 \\
       k_{31}^0 & k_{32}^0 & k_{33}^0 & k_{34}^0 \\
       k_{41}^0 & k_{42}^0 & k_{43}^0 & k_{44}^0
     \end{bmatrix}}\begin{matrix}0\\1\\2\\3\end{matrix} \longrightarrow
     \stackrel{\begin{matrix}
         \hphantom{k_{11}^0} &
         \hphantom{k_{12}^0} &
         \hphantom{k_{13}^0} &
         \hphantom{k_{14}^0} &
         \hphantom{0} &
         \hphantom{0} \\ 0 & 1 & 2 & 3 & 4 & 5 \end{matrix}}{
     \begin{bmatrix}
       k_{11}^0 & k_{12}^0 & k_{13}^0 & k_{14}^0 & 0 & 0 \\
       k_{21}^0 & k_{22}^0 & k_{23}^0 & k_{24}^0 & 0 & 0\\
       k_{31}^0 & k_{32}^0 & k_{33}^0 & k_{34}^0 & 0 & 0\\
       k_{41}^0 & k_{42}^0 & k_{43}^0 & k_{44}^0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0\\
       0 & 0 & 0 & 0 & 0 & 0
     \end{bmatrix}}\begin{matrix}0\\1\\2\\3\\4\\5\\ \end{matrix}

The local-to-global-DOF mapping can be generalized as:

.. math::

   \text{GLOBAL DOF} = \text{NUMBER OF DOF PER NODE}\times\text{NODE
       NUMBER} + \text{LOCAL DOF}

:math:`\text{LOCAL DOF}` is one of the following integers:

+--------------------+-------------+
| DOF description    | LOCAL DOF   |
+====================+=============+
| :math:`x`          | 0           |
+--------------------+-------------+
| :math:`y`          | 1           |
+--------------------+-------------+
| :math:`z`          | 2           |
+--------------------+-------------+
| :math:`\theta_x`   | 3           |
+--------------------+-------------+
| :math:`\theta_y`   | 4           |
+--------------------+-------------+
| :math:`\theta_z`   | 5           |
+--------------------+-------------+
| :math:`T`          | 6           |
+--------------------+-------------+

.. warning::

   This local-to-global-DOF mapping is only applicable for 0 based array
   indexing. For 1 based indexing, the mapping must be adjusted.

The element freedom table for element 1 is ``eft=[0,1,4,5]`` and
:math:`\left[k^1\right]` and :math:`\left[K\right]` are:

.. math::

     \stackrel{\begin{matrix}
         \hphantom{k_{11}^0} &
         \hphantom{k_{12}^0} &
         \hphantom{k_{13}^0} &
         \hphantom{k_{14}^0} \\ 0 & 1 & 4 & 5 \end{matrix}}{
     \begin{bmatrix}
       k_{11}^1 & k_{12}^1 & k_{13}^1 & k_{14}^1 \\
       k_{21}^1 & k_{22}^1 & k_{23}^1 & k_{24}^1 \\
       k_{31}^1 & k_{32}^1 & k_{33}^1 & k_{34}^1 \\
       k_{41}^1 & k_{42}^1 & k_{43}^1 & k_{44}^1
     \end{bmatrix}}\begin{matrix}0\\1\\4\\5\end{matrix} \longrightarrow
     \stackrel{\begin{matrix}
         \hphantom{k_{11}^0+k_{11}^1} &
         \hphantom{k_{12}^0+k_{11}^1} &
         \hphantom{k_{13}^0} &
         \hphantom{k_{14}^0} &
         \hphantom{k_{11}^1} &
         \hphantom{k_{11}^1} \\ 0 & 1 & 2 & 3 & 4 & 5 \end{matrix}}{
     \begin{bmatrix}
       k_{11}^0+k_{11}^1 & k_{12}^0+k_{11}^1 & k_{13}^0 & k_{14}^0 & k_{13}^1 & k_{23}^1 \\
       k_{21}^0+k_{12}^1 & k_{22}^0+k_{22}^1 & k_{23}^0 & k_{24}^0 & k_{14}^1 & k_{24}^1 \\
       k_{31}^0 & k_{32}^0 & k_{33}^0 & k_{34}^0 & 0 & 0\\
       k_{41}^0 & k_{42}^0 & k_{43}^0 & k_{44}^0 & 0 & 0\\
       k_{31}^1 & k_{32}^1 & 0 & 0 & k_{33}^1 & k_{34}^1\\
       k_{41}^1 & k_{42}^1 & 0 & 0 & k_{43}^1 & k_{44}^1\\
     \end{bmatrix}}\begin{matrix}0\\1\\2\\3\\4\\5\\ \end{matrix}

.. note::

   When assembling the global finite element stiffness, the
   following questions must be answered in the affirmitive:

   -  Is the value of :math:`K_{II}` equal to the sum of all elements
      connected to node :math:`I`?

   -  Is the value of :math:`K_{IJ}`, :math:`I\ne J`, equal to the negative
      of all elements connecting nodes :math:`I` and :math:`J`?

   -  Is the stiffness symmetric?

Looking at the global system in , it is clear that each of the three
preceding questions can be answered in the affirmitive.

Assembler Computational Implementation
--------------------------------------

Algorithms for assemblers for the global stiffness and global force are presented.

Global Stiffness
~~~~~~~~~~~~~~~~

The algorithm for assembling the global stiffness is

#. **TASK**: Assemble the global stiffness ``K``

   #. **INITIALIZE**

      - **SET** ``N`` = number of nodes in the mesh
      - **SET** ``m`` = number of degrees of freedom per node
      - **SET** ``K[N*m,N*m]`` = 0

   #. **FOR EACH** element **DO**

      - **SET** ``Ke`` = element stiffness
      - **SET** ``n`` = number of nodes on element
      - **SET** ``nodes[n]`` = node numbers of the element
      - **SET** ``eft[n*m]`` = 0
      - **SET** ``i`` = 0

      **COMMENT**: Construct the element freedom table

      #. **DO WHILE** ``i`` < ``n``

         - **SET** ``ni`` = ``nodes[i]``
	 - **SET** ``j`` = ``0``

	 #. **DO WHILE** ``j`` < ``m``:

	    - **SET** ``eft[i+j]`` = ``ni*m`` + ``j``
	    - **SET** ``j`` = ``j`` + 1

	 - **SET** ``i`` = ``i`` + 1

      **COMMENT** Map the local degrees of freedom to the global

      - **SET** ``i`` = 0

      #. **DO WHILE** ``i`` < ``n*m``

	 - **SET** ``I`` = ``eft[i]``
         - **SET** ``j`` = ``i``

	 #. **DO WHILE** ``j`` < ``n*m``

	    - **SET** ``J`` = ``eft[j]``
            - **SET** ``K[I,J]`` = ``K[I,J]`` + ``Ke[i,j]``
            - **SET** ``K[J, I]`` = ``K[I, J]``
	    - **SET** ``j`` = ``j`` + 1

         - **SET** ``i`` = ``i`` + 1

In ``pyfem2``, assembly of the global stiffness is performed in the method
``FiniteElementModel.assemble_global_stiffness_udof`` (``udof`` stands for
"uniform degree of freedom").

The heart of the assembly process is constructing the element freedom table,
represented by the symbol ``eft``. The element freedom table contains the global
degree of freedom number for the nodes in the element. For elements having
uniform degrees of freedom, this table can be constructed on the fly using, as shown.

Stiffness assembler example
...........................

A representative python code that performs stiffness assembly is the function
``AssemblePlaneTrussStiffness``, shown below. ``AssemblePlaneTrussStiffness`` is invoked as

.. code:: python

   K = AssemblePlaneTrussStiffness(coord, elecon, elemat, elefab)

The arguments ``AssemblePlaneTrussStiffness`` to are:

+--------------+------------------------------------------------------------------------------------------------------+
| ``coord``    | Nodal coordinates. ``coord[n,i]`` is the ith coordinate of node n.                                   |
+--------------+------------------------------------------------------------------------------------------------------+
| ``elecon``   | Element connectivity. ``elecon[e,n]`` is the nth internal node label of element e.                   |
+--------------+------------------------------------------------------------------------------------------------------+
| ``elemat``   | Element material properties. ``elemat[e]`` are the material properties of element e.                 |
+--------------+------------------------------------------------------------------------------------------------------+
| ``elefab``   | Element fabrication properties. ``elefab[e]`` are the element fabrication properties of element e.   |
+--------------+------------------------------------------------------------------------------------------------------+

The output is:

+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``K``   | Global stiffness matrix stored as a full (N\*m, N\*m) symmetric matrix, where m is the number of degrees of freedom per node and N the total number of nodes. |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. code:: python

    def AssemblePlaneTrussStiffness(coord, elecon, elemat, elefab):
	coord = asarray(coord)
	elecon = asarray(elecon)

        # number of nodes and degrees of freedom per node
        m = 2
        N = coord.shape[0]
        numele = elecon.shape[0]

        # compute the element stiffness and scatter to global array
        K = zeros((N*m, N*m))
        for e in range(numele):

	    nodes = elecon[e]  # element connectivity
	    x = coord[nodes]  # nodal coordinates
	    n = len(nodes)  # number of nodes on element

            # Element stiffness
            Ke = PlaneTrussElementStiffness(x, elemat[e], elefab[e])

            # Element freedom table
	    eft = zeros(n*m, dtype=int)
	    for i in range(n):
	        ni = nodes[i]
		for j in range(m):
		    eft[i+j] = ni*m + j

            # Merge element stiffness in to global
            for i in range(n*m):
                I = eft[i]
                for j in range(i, n*m):
                    J = eft[j]
                    K[I, J] += Ke[i,j]
                    K[J, I] = K[I, J]
        return K

The subordinate function ``PlaneTrussElementStiffness`` computes the element stiffness matrix (see :ref:`BarElements`):

.. code:: python

    def PlaneTrussElementStiffness(x, E, A):
        v = x[1] - x[0]
        h = sqrt(dot(v, v))
        n = v / h
        nn = outer(n, n)

        # Assemble element stiffness
        i, j = 2, 4
        k = zeros((4, 4))
        k[0:i, 0:i] = k[i:j, i:j] =  nn # upper left and lower right 2x2
        k[0:i, i:j] = k[i:j, 0:i] = -nn # lower left and upper right 2x2
        return A * E / h * k

.. note::

   ``AssemblePlaneTrussStiffness`` is valid only for plane truss elements.

.. note::

  Clarity of code was chosen over performance in the implementation of
  ``AssemblePlaneTrussStiffness``. The process of creating the element freedom
  table and merging element equations in to the global equations can be
  optimized greatly.

Global Force
~~~~~~~~~~~~

The algorithm for assembling the global force follows closely that of assembling the global stiffness:

#. **TASK**: Assemble the global force ``F``

   #. **INITIALIZE**

      - **SET** ``N`` = number of nodes in the mesh
      - **SET** ``m`` = number of degrees of freedom per node
      - **SET** ``F[N*m]`` = 0

   #. **FOR EACH** element **DO**

      - **SET** ``Fe`` = element force
      - **SET** ``n`` = number of nodes on element
      - **SET** ``nodes[n]`` = node numbers of the element
      - **SET** ``eft[n*m]`` = 0
      - **SET** ``i`` = 0

      **COMMENT**: Construct the element freedom table

      #. **DO WHILE** ``i`` < ``n``

         - **SET** ``ni`` = ``nodes[i]``
	 - **SET** ``j`` = ``0``

	 #. **DO WHILE** ``j`` < ``m``:

	    - **SET** ``eft[i+j]`` = ``ni*m`` + ``j``
	    - **SET** ``j`` = ``j`` + 1

	 - **SET** ``i`` = ``i`` + 1

      **COMMENT** Map the local degrees of freedom to the global

      - **SET** ``i`` = 0

      #. **DO WHILE** ``i`` < ``n*m``

	 - **SET** ``I`` = ``eft[i]``
         - **SET** ``F[I]`` = ``F[I]`` + ``Fe[i]``

In ``pyfem2``, assembly of the global stiffness is performed in the method
``FiniteElementModel.assemble_global_force_udof``.
