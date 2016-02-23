.. _TrussModel:

TrussSolution: A Complete Truss Program
=======================================

-  Chapter [sec:first],

-  Chapter [sec:intro.conventions],

-  Chapter [ch:bc\_and\_load],

-  Chapter [sec:ass],

-  Chapter [sec:appbc],

-  `IFEM Chapter 2, “The Direct Stiffness Method I” <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch02.d/IFEM.Ch02.pdf>`__
-  `IFEM Chapter 3, “The Direct Stiffness Method II” <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch03.d/IFEM.Ch03.pdf>`__

Overview
--------

models truss-like structures, like the one shown in Figure
[fig:truss.intro] in one, two, and three-dimensions.

-  models truss members as two-node elastic bars having constant
   cross-sectional area and Young’s modulus,

-  supports truss members defined in 1, 2, or 3 dimensions,

-  uses the direct stiffness method to assemble and solve the global
   system of equations,

-  supports prescribed displacements and nodal forces,

-  uses linear solver to solve the global system of equations, and

-  writes results to `vtk <www.vtk.ort>`__ formatted files.

.. figure:: IntroTruss.png
   :align: center

   Finite element analysis of a truss-like structure

Program limitations
~~~~~~~~~~~~~~~~~~~

-  is valid only for truss-like structures whose members are two-node
   elastic bars,

-  does not support distributed loads,

-  does not support multi-freedom constraints, and

-  sacrifices speed and efficiency for clarity.

Element types
-------------

supports the following element types:

-  : one-dimensional two-node elastic bar.

-  : two-dimensional two-node elastic bar.

-  : three-dimensional two-node elastic bar.

The different element types are described in Appendix [app:elem.lib].

Truss program design
--------------------

A callgraph of the truss program, shown in Figure [fig:truss.callgraph],
demonstrates its overal design. As described in “Chapter 1: Introduction
and Conventions”, the program is split broadly in to preprocessing,
processing, and postprocessing stages. The user defines the model
definition in the preprocessing stage and the truss program performs the
processing and postprocessing steps.

.. figure:: TrussCallGraph.png
   :align: center

   Callgraph of Truss FE program

The truss program script
~~~~~~~~~~~~~~~~~~~~~~~~

| The function , shown in Listing [lst:truss.app], is an executive type
  function that organizes the processing and postprocessing analysis
  stages. It is invoked as

.. code:: python

    u, R, p, s = TrussSolution(nodtab, eletab, elemat, elefab, bcs, cloads)

The arguments to are

+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``nodtab``   | Table of defining nodes. is the node label of node and are the coordinates of node .                                                                                                                                                                                    |
+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``eletab``   | Table of defining elements. is the element label of element , is the element type, and are the node labels of nodes forming element .                                                                                                                                   |
+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``elemat``   | Element material data. For the Link2 element, the only element material property is the Young’s modulus E. is given as either a float if all elements have the same constant Young’s modulus or a list , where is the total number of elements.                         |
+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``elefab``   | Element fabrication data. For the Link2 element, the only element fabrication property is the cross-sectional area A. is given as either a float if all elements have the same constant cross sectional area or a list , where is the total number of elements.         |
+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``bcs``      | Displacement boundary condition specification. List of lists defining the prescribed displacements in the model. Each sublist contains three elements defining the node ID[s] constrained, the degree[s] of freedom constrained, and the magnitude of the constraint.   |
+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``cloads``   | Concentrated (nodal) force specification. The composition of is identical to that of .                                                                                                                                                                                  |
+--------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The outputs from are

+---------+------------------------------------------------------------------+
| ``u``   | Nodal displacements. , where is the total number of nodes.       |
+---------+------------------------------------------------------------------+
| ``R``   | Nodal reactions. , where is the total number of nodes.           |
+---------+------------------------------------------------------------------+
| ``p``   | Element axial forces. , where is the total number of elements.   |
+---------+------------------------------------------------------------------+
| ``s``   | Element stresses , where is the total number of elements.        |
+---------+------------------------------------------------------------------+

.. code:: python

    def TrussSolution(nodtab,eletab,elemat,elefab,bcs,cload,jobid='Job-1'):
        # Format user input
        nodmap,coord,elemap,eletyp,elecon=FormatNodesAndElements(nodtab,eletab)
        # Truss program valid only for Link elements
        if len(set(eletyp)) != 1:
            raise ValueError('TrussSolution expected only one element type')
        if any([et not in (L1D2, L2D2, L3D2) for et in eletyp]):
            raise ValueError('TrussSolution expected element types '
                             'L1D2, L2D2, or L3D2')
        # number of elements, nodes
        numele, numnod = elecon.shape[0], coord.shape[0]
        # number of degrees of freedom per node
        ndof = 1 if coord.ndim == 1 else coord.shape[1]
        # make sure properties are arrays
        if not IsListlike(elemat): elemat = array([elemat] * numele)
        if not IsListlike(elefab): elefab = array([elefab] * numele)
        # Format the boundary conditions
        doftags, dofvals = FormatBoundaryConditions(coord, nodmap, eletyp, elecon,
                                                    bcs, cload)
        # Limit doftags, dofvals to just the DOF of interest
        doftags, dofvals = doftags[:, :ndof], dofvals[:, :ndof]
        # Assemble the global stiffness and force
        K = AssembleGlobalStiffnessUDOF(coord, elemap, eletyp, elecon, elemat, elefab)
        F, Q = AssembleGlobalForceUDOF(coord, elemat, eletyp, elecon, doftags, dofvals)
        Kbc, Fbc = ApplyBoundaryConditionsUDOF(K, F+Q, doftags, dofvals)
        try:
            u = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')
        # Total force, including reaction, and reaction
        Ft = dot(K, u)
        R = Ft - F - Q
        # reshape u, R to be the same shape as coord
        u, R = u.reshape(coord.shape), R.reshape(coord.shape)
        p = TrussIntForces(coord, elecon, elemat, elefab, u)
        s = TrussStresses(p, elefab)
        if jobid is not None:
            WriteFEResults(jobid, coord, nodmap, elemap, eletyp, elecon,
                           u=u, p=p, R=R, s=s)
        return u, R, p, s

Example model driver
^^^^^^^^^^^^^^^^^^^^

Consider the space truss shown in Figure [fig:truss.space\_truss], the
model driver for this problem is given Listing [lst:truss.model].

.. figure:: SpaceTruss1.jpeg
   :align: center

   Example space truss

.. code:: python

    nodtab = [[1, 72, 0, 0], [2, 0, 36, 0], [3, 0, 36, 72], [4, 0, 0, -48]]
    eletab = [[1, 1, 2], [2, 1, 3], [3, 1, 4]]
    A = [.302, .729, .187]
    E = 10e6

    # Boundary conditions
    bcs = [[1, Y, 0]], [[2,3,4], (X,Y,Z), 0]]

    # Concentrated force in 'z' direction on node 1
    cloads = [[1, Z, -1000]]

    u, R, p, s = TrussSolution(nodtab, eletab, E, A, bcs, cloads)

Displacement solution and recovery of reactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| The linear system
  :math:`\boldsymbol{K}^*\boldsymbol{u}=\boldsymbol{F}^*`, where
  :math:`\boldsymbol{K}^*` and :math:`\boldsymbol{F}^*` are the boundary
  condition modified global stiffness and force, respectively, is solved
  using the function that is a part of the module. is invoked as

.. code:: python

      u = numpy.linalg.solve(Kbc, Fbc)

Reaction forces :math:`\boldsymbol{R}` are recovered as

.. math:: \boldsymbol{R} = \boldsymbol{K}\boldsymbol{u} - \boldsymbol{F}

where :math:`\boldsymbol{K}` and :math:`\boldsymbol{F}` are the global
stiffness and force, respectively (before application of boundary
conditions), and :math:`\boldsymbol{u}` is the displacement found from
the call to .

Postprocessing
~~~~~~~~~~~~~~

| With nodal displacements known, postprocessing can begin. In this
  program, postprocessing consists of determining the element internal
  forces and stresses. computes the axial internal forces of truss
  members and is invoked as

.. code:: python

      p = TrussIntForces(coor, nodmap, elecon, elemat, elefab, u)

The arguments , , , , , and are described in earlier sections. uses the
subordinate function to determine the internal force in individual
elements.

| The internal stress is computed in the function and is invoked as

.. code:: python

    s = TrussStresses(p, elefab)

The ouput from is

+---------+-------------------------------------------------------------------------------+
| ``p``   | Element axial forces , where is the total number of elements. Computed by .   |
+---------+-------------------------------------------------------------------------------+

.. code:: python

    def TrussIntForces(coor, elecon, elemat, elefab, u):
        p = zeros(elecon.shape[0])
        for (e, c) in enumerate(elecon):
            c = [nodmap[n] for n in c]
            p[e] = Link2IntForce(coor[c], elemat[e], elefab[e], u[c])
        return p

    def Link2IntForce(xc, E, A, uc):
        x = xc[1] - xc[0]
        u = uc[1] - uc[0]
        Xu = dot(x, u)
        L = sqrt(dot(x, x))
        return E * A / L * Xu / L

.. code:: python

    def TrussStresses(p, elefab):
        return p / elefab

Writing FE results
~~~~~~~~~~~~~~~~~~

writes the results to a `vtk <www.vtk.org>`__ finite element database.
files and are viewable in several commercial and open source
visualization products, including the open source
`ParaView <http://www.paraview.org>`__. is invoked as

.. code:: python

    WriteFEResults(jobid, coord, nodmap, elemap, eletyp, elecon, u=u, **kwds)

The arguments to are

+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``jobid``    | String identifying the simulation                                                                                                        |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``coord``    | Nodal coordinates , where is the total number of nodes. Node IDs are implied by the row number.                                          |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``elemap``   | Mapping from external node label to internal node number. .                                                                              |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``eletyp``   | Element type. is the element type of element e.                                                                                          |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``elecon``   | Element connectivity (nodes defining ends of element). is the th node ID of the th element. Element IDs are implied by the row number.   |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``u``        | Nodal displacements , where is the total number of nodes. Node IDs are implied by the row number.                                        |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+
| ``**kwds``   | Keyword arguments. Pass nodal and element data as keywords.                                                                              |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------+

The argument to allows one to store arbitrary nodal and element results
to the output database. For example, to store the nodal reactions,
element forces, and element stresses, would be invoked as

Examples
--------

Example 1
~~~~~~~~~

Consider the bridge shown in Figure [fig.truss.bridge]. The user model
definition script is

.. figure:: Bridge.png
   :align: center

   Six-bay bridge plane truss. (a) physical problem; (b) finite element idealization. (Credit: `Ch. 21 <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch21.d/IFEM.Ch21.pdf>`__ of Prof. Felippa’s introductory finite element course materials)

.. code:: python

    coor = [[1,0,0,0], [2,10,5,0], [3,10,0,0], [4,20,8,0], [5,20,0,0],
            [6,30,9,0], [7,30,0,0], [8,40,8,0], [9,40,0,0], [10,50,5,0],
            [11,50,0,0],[12,60,0,0]]
    elecon = [[1,1,3], [2,3,5], [3,5,7], [4,7,9], [5,9,11], [6,11,12],
              [7,1,2], [8,2,4], [9,4,6], [10,6,8], [11,8,10], [12,10,12],
              [13,2,3], [14,4,5], [15,6,7], [16,8,9], [17,10,11], [18,2,5],
              [19,4,7], [20,7,8], [21,9,10]]
    E = 1000
    Abot, Atop, Abat, Adia = 2, 10, 3, 1
    A = [Abot, Abot, Abot, Abot, Abot, Abot,
         Atop, Atop, Atop, Atop, Atop, Atop,
         Abat, Abat, Abat, Abat, Abat,
         Adia, Adia, Adia, Adia]
    cloads = [[3,Y,-10.],[5,Y,-10.],[7,Y,-16.],[9,Y,-10.],[11,Y,-10.]]
    bcs = [[1, ALL, 0.], [12, Y, 0.], [ALL, Z, 0]]
    u, R, p, s = TrussSolution(nodtab, eletab, E, A, bcs, cloads)

The computed displacement and nodal reactions are

.. code:: python

    nodal displacements               nodal reactions
    [[ 0.       0.       0.     ]     [[  0.  28.   0.]
     [ 0.80954 -1.7756   0.     ]      [  0.   0.   0.]
     [ 0.28    -1.79226  0.     ]      [  0.   0.   0.]
     [ 0.899   -2.29193  0.     ]      [  0.   0.   0.]
     [ 0.56    -2.3166   0.     ]      [  0.   0.   0.]
     [ 0.8475  -2.38594  0.     ]      [ -0.  -0.   0.]
     [ 0.8475  -2.42194  0.     ]      [ -0.   0.   0.]
     [ 0.796   -2.29193  0.     ]      [ -0.   0.   0.]
     [ 1.135   -2.3166   0.     ]      [  0.  -0.   0.]
     [ 0.88546 -1.7756   0.     ]      [  0.   0.   0.]
     [ 1.415   -1.79226  0.     ]      [  0.   0.   0.]
     [ 1.695    0.       0.     ]]     [  0.  28.   0.]]

The undeformed and deformed plots, generated by ParaView, are shown in
Figure [fig.truss.bridge\_res]

.. figure:: TrussExample1.png
   :align: center

   Undeformed and deformed plots of the truss shown in Figure
   [fig:truss.bridge]. The deformed plots show contours of :math:`y`
   displacement and magnitude of the reaction forces.

Example 2
~~~~~~~~~

.. figure:: TrussExample2.jpg
   :align: center

   Truss example 2

.. code:: python

    nodtab = [[1, 0, 0], [2, 3, 4], [3, 0, 4]]
    eletab = [[1, 1, 2], [2, 1, 3]])
    E, A = 70e9, 5 * .01 * .01
    bcs = [[1, X, -.05], [[2,3], (X,Y), 0]]
    cloads = [[1, Y, 1000e3]]
    u, R, p, s = TrussSolution(nodtab, eletab, E, A, bcs, cloads)

The nodal displacements are

.. code:: python

    [[-0.05     0.08828]
     [ 0.       0.     ]
     [ 0.       0.     ]]

Example 3
~~~~~~~~~

Redo Example 2, but in 3D.

.. code:: python

    nodtab = [[1,0,0,0], [2,3,4,0], [3,0,4,0]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    E, A = 70e9, 5 * .01 * .01
    bcs = [[1, X, -.05], [1, Z, 0], [[2,3], (X,Y,Z), 0]]
    cloads = [[1, Y, 1000e3]]
    u, R, p, s = TrussSolution(nodtab, eletab, E, A, bcs, cloads)

The nodal displacements are

.. code:: python

    [[-0.05     0.08828  0.     ]
     [ 0.       0.       0.     ]
     [ 0.       0.       0.     ]]
