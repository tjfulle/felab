
.. _Truss2:

Space truss
===========

Objective
---------

This example demonstrates the use of the ``TrussModel``  to solve for the deflections in a space truss.  ``TrussModel`` is described in more detail in :ref:`TrussModel`.

Problem description
-------------------

Consider the six bay plane truss shown (credit: `Ch. 21 <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch21.d/IFEM.Ch21.pdf>`__ of Prof. Felippa’s introductory finite element course materials).

.. figure:: SpaceTrussEx.png
   :align: center

   25-member space truss model of a transmission tower. (a): Geometry definition and node numbers; (b) element numbers.

Model script
------------

The problem is defined and solved in ``pyfem2`` as follows:

.. code:: python

   """
   pyfem2 tutorial demo program: Plane truss.

   """

   from numpy import *
   from pyfem2 import *

   V = TrussModel()

   # Create mesh and define function space
   nodtab = [[1,-37.5,0,200],[2,37.5,0,200],[3,-37.5,37.5,100],
             [4,37.5,37.5,100],[5,37.5,-37.5,100],[6,-37.5,-37.5,100],
             [7,-100,100,0],[8,100,100,0],[9,100,-100,0],[10,-100,-100,0]]
   eletab = [[ 1, 1, 2],[ 2, 1, 4],[ 3, 2, 3],[ 4, 1, 5],[ 5, 2, 6],
             [ 6, 2, 4],[ 7, 2, 5],[ 8, 1, 3],[ 9, 1, 6],[10, 3, 6],
             [11, 4, 5],[12, 3, 4],[13, 5, 6],[14, 3,10],[15, 6, 7],
             [16, 4, 9],[17, 5, 8],[18, 4, 7],[19, 3, 8],[20, 5,10],
             [21, 6, 9],[22, 6,10],[23, 3, 7],[24, 5, 9],[25, 4, 8]]
   V.Mesh(nodtab=nodtab, eletab=eletab)

   # Define element blocks
   V.ElementBlock('ElemenbBlock1', ALL)

   # Assign properties
   V.Material('Material-1')
   V.materials['Material-1'].Elastic(E=10e6, Nu=.333)
   A = [0.033, 2.015, 2.015, 2.015, 2.015, 2.823, 2.823, 2.823, 2.823, 0.01,
        0.01, 0.014, 0.014, 0.98, 0.98, 0.98, 0.98, 1.76, 1.76, 1.76, 1.76,
        2.44, 2.44, 2.44, 2.44]
   V.AssignProperties('ElemenbBlock1', Link3D2, 'Material-1', A=A)

   # Define boundary conditons
   V.FixNodes((7, 8, 9, 10))

   # Define concentrated loads
   P1, P2, P3, P4 = 1000, 10000, -5000, 500
   V.ConcentratedLoad(1, X, P1)
   V.ConcentratedLoad(1, Y, P2)
   V.ConcentratedLoad(1, Z, P3)
   V.ConcentratedLoad(2, Y, P2)
   V.ConcentratedLoad(2, Z, P3)
   V.ConcentratedLoad((3, 6), X, P4)

   # Solve and write results
   V.Solve()
   V.WriteResults('Truss2.exo')

How does it work?
-----------------

The complete code can be found in the files ``tutorials/Truss2.py``. We now examine the preceding program in detail.

The first lines of the program,

.. code:: python

   from numpy import *
   from pyfem2 import *

import objects from the ``numpy`` and ``pyfem2`` namespaces in to the program.
`numpy <http://www.numpy.org>`__ is a python package providing numerical data
types and procedures. The key imports from the ``pyfem2`` library is
the ``TrussModel``.

The statement

.. code:: python

   V = TrussModel()

creates the finite element model.  The finite element mesh is created by defining tables of nodes and elements (see :ref:`NodeDefinition` and :ref:`ElementDefinition`) and passing them to the ``Mesh`` method:

.. code:: python

   nodtab = [[1,-37.5,0,200],[2,37.5,0,200],[3,-37.5,37.5,100],
             [4,37.5,37.5,100],[5,37.5,-37.5,100],[6,-37.5,-37.5,100],
             [7,-100,100,0],[8,100,100,0],[9,100,-100,0],[10,-100,-100,0]]
   eletab = [[ 1, 1, 2],[ 2, 1, 4],[ 3, 2, 3],[ 4, 1, 5],[ 5, 2, 6],
             [ 6, 2, 4],[ 7, 2, 5],[ 8, 1, 3],[ 9, 1, 6],[10, 3, 6],
             [11, 4, 5],[12, 3, 4],[13, 5, 6],[14, 3,10],[15, 6, 7],
             [16, 4, 9],[17, 5, 8],[18, 4, 7],[19, 3, 8],[20, 5,10],
             [21, 6, 9],[22, 6,10],[23, 3, 7],[24, 5, 9],[25, 4, 8]]
   V.Mesh(nodtab=nodtab, eletab=eletab)

The finite element model ``V`` requires that every element in the mesh be
assigned an element type and material constitutive relations. The assignment
occurs by grouping elements of the same type together in to element blocks and
then assigning to those element blocks material and fabrication properties.
For this problem, every element in the mesh is a ``Link3D2`` two-node elastic
bar element (an elastic bar element with two degrees of freedom per node). The
element block, named ``ElementBlock1``, containing all of the elements in the
mesh is created by:

.. code:: python

   V.ElementBlock('ElemenbBlock1', ALL)

The symbol ``ALL`` informs the ``ElementBlock`` method to assign all elements
in the mesh to the block ``ElementBlock1``. Material properties are defined by
the ``Material`` method and elements in a block are assigned material and
fabrication properties by the ``AssignProperties`` method:

.. code:: python

   V.Material('Material-1')
   V.materials['Material-1'].Elastic(E=10e6, Nu=.333)
   A = [0.033, 2.015, 2.015, 2.015, 2.015, 2.823, 2.823, 2.823, 2.823, 0.01,
        0.01, 0.014, 0.014, 0.98, 0.98, 0.98, 0.98, 1.76, 1.76, 1.76, 1.76,
        2.44, 2.44, 2.44, 2.44]
   V.AssignProperties('ElemenbBlock1', Link3D2, 'Material-1', A=A)

The method ``AssignProperties`` takes as input the name of the element block
to which properties are being assigned, the element type for elements in the
block, the material model name, and any element fabrication properties. For
``Link2D2`` elements, the area ``A`` is the only fabrication property.

The next step is to specify the boundary conditions :math:`u_{7}=u_{8}=u_9=u_{10}=0`:

.. code:: python

   V.FixNodes((7, 8, 9, 10))

The point forces at nodes 1, 2, 3, and 6

.. code:: python

   P1, P2, P3, P4 = 1000, 10000, -5000, 500
   V.ConcentratedLoad(1, X, P1)
   V.ConcentratedLoad(1, Y, P2)
   V.ConcentratedLoad(1, Z, P3)
   V.ConcentratedLoad(2, Y, P2)
   V.ConcentratedLoad(2, Z, P3)
   V.ConcentratedLoad((3, 6), X, P4)

Finally, the unknown displacements are determined by solving the model and the model results are written to an ExodusII output file

   V.Solve()
   V.WriteResults('Truss2.exo')

The deformed geometry, viewed in `ParaView <http://www.paraview.org>`__, is shown below

.. figure:: SpaceTrussEx2.png
   :align: center

The deformed plots show contours of displacement magnitude.
