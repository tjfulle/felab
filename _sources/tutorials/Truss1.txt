
.. _Truss1:

Plane truss
===========

Objective
---------

This example demonstrates the use of the ``StaticStep`` method of the ``FiniteElementModel``  to solve for the deflections in a plane truss.

Problem description
-------------------

Consider the six bay plane truss shown (credit: `Ch. 21 <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch21.d/IFEM.Ch21.pdf>`__ of Prof. Felippa’s introductory finite element course materials).

.. figure:: Bridge.png
   :align: center

   Six-bay bridge plane truss. (a) physical problem; (b) finite element idealization.

Model script
------------

The problem is defined and solved in ``pyfem2`` as follows:

.. code:: python

   """
   pyfem2 tutorial demo program: Plane truss.

   """

   from numpy import *
   from pyfem2 import *

   # Create the model problem
   V = FiniteElementModel(jobid='Truss1')

   # Create the mesh from tables of nodes and elements
   nodtab = [[1,0,0], [2,10,5], [3,10,0], [4,20,8], [5,20,0],
             [6,30,9], [7,30,0], [8,40,8], [9,40,0], [10,50,5],
             [11,50,0],[12,60,0]]
   eletab = [[1,1,3], [2,3,5], [3,5,7], [4,7,9], [5,9,11], [6,11,12],
             [7,1,2], [8,2,4], [9,4,6], [10,6,8], [11,8,10], [12,10,12],
             [13,2,3], [14,4,5], [15,6,7], [16,8,9], [17,10,11], [18,2,5],
             [19,4,7], [20,7,8], [21,9,10]]
   V.Mesh(nodtab=nodtab, eletab=eletab)

   # Define an element block of 3D 2-node link elements
   V.ElementBlock('ElementBlock1', ALL)

   # Create a material and define the elastic properties
   mat = V.Material('Material-1', elastic={'E':1000, 'Nu':.29})
   Abot, Atop, Abat, Adia = 2, 10, 3, 1
   A = [Abot, Abot, Abot, Abot, Abot, Abot,
        Atop, Atop, Atop, Atop, Atop, Atop,
        Abat, Abat, Abat, Abat, Abat,
        Adia, Adia, Adia, Adia]
   V.AssignProperties('ElementBlock1', ElasticLink2D2, mat, A=A)

   # Define a load step
   step = V.StaticStep()

   # Apply boundary conditions
   step.PrescribedBC(1, (X,Y))
   step.PrescribedBC(12, Y)

   # Apply concentrated loads
   step.ConcentratedLoad((3,5,9,11), Y, -10)
   step.ConcentratedLoad(7, Y, -16)

   step.run()

The computed displacement and nodal reactions are

.. code:: python

    nodal displacements               nodal reactions
    [[ 0.       0.       ]            [[  0.  28.]
     [ 0.80954 -1.7756   ]             [  0.   0.]
     [ 0.28    -1.79226  ]             [  0.   0.]
     [ 0.899   -2.29193  ]             [  0.   0.]
     [ 0.56    -2.3166   ]             [  0.   0.]
     [ 0.8475  -2.38594  ]             [ -0.  -0.]
     [ 0.8475  -2.42194  ]             [ -0.   0.]
     [ 0.796   -2.29193  ]             [ -0.   0.]
     [ 1.135   -2.3166   ]             [  0.  -0.]
     [ 0.88546 -1.7756   ]             [  0.   0.]
     [ 1.415   -1.79226  ]             [  0.   0.]
     [ 1.695    0.       ]]            [  0.  28.]]


How does it work?
-----------------

The complete code can be found in the files ``pyfem2/data/Truss1.py``. We now examine the preceding program in detail.

The first lines of the program,

.. code:: python

   from numpy import *
   from pyfem2 import *

import objects from the ``numpy`` and ``pyfem2`` namespaces in to the program.
`numpy <http://www.numpy.org>`__ is a python package providing numerical data
types and procedures. The key imports from the ``pyfem2`` library is
the ``FiniteElementModel``.

The statement

.. code:: python

   V = FiniteElementModel(jobid='Truss1')

creates the finite element model.  The finite element mesh is created by defining tables of nodes and elements (see :ref:`NodeDefinition` and :ref:`ElementDefinition`) and passing them to the ``Mesh`` method:

.. code:: python

   nodtab = [[1,0,0], [2,10,5], [3,10,0], [4,20,8], [5,20,0],
             [6,30,9], [7,30,0], [8,40,8], [9,40,0], [10,50,5],
             [11,50,0],[12,60,0]]
   eletab = [[1,1,3], [2,3,5], [3,5,7], [4,7,9], [5,9,11], [6,11,12],
             [7,1,2], [8,2,4], [9,4,6], [10,6,8], [11,8,10], [12,10,12],
             [13,2,3], [14,4,5], [15,6,7], [16,8,9], [17,10,11], [18,2,5],
             [19,4,7], [20,7,8], [21,9,10]]
   V.Mesh(nodtab=nodtab, eletab=eletab)

The finite element model ``V`` requires that every element in the mesh be
assigned an element type and material constitutive relations. The assignment
occurs by grouping elements of the same type together in to element blocks and
then assigning to those element blocks material and fabrication properties.
For this problem, every element in the mesh is a ``ElasticLink3D2`` two-node
elastic bar element (an elastic bar element with two degrees of freedom per
node). The element block, named ``ElementBlock1``, containing all of the
elements in the mesh is created by:

.. code:: python

   V.ElementBlock('ElementBlock1', ALL)

The symbol ``ALL`` informs the ``ElementBlock`` method to assign all elements
in the mesh to the block ``ElementBlock1``. Material properties are defined by
the ``Material`` method and elements in a block are assigned material and
fabrication properties by the ``AssignProperties`` method:

.. code:: python

   mat = V.Material('Material-1')
   mat.Elastic(E=1000, Nu=.29)
   Abot, Atop, Abat, Adia = 2, 10, 3, 1
   A = [Abot, Abot, Abot, Abot, Abot, Abot,
        Atop, Atop, Atop, Atop, Atop, Atop,
        Abat, Abat, Abat, Abat, Abat,
        Adia, Adia, Adia, Adia]
   V.AssignProperties('ElementBlock1', ElasticLink2D2, mat, A=A)

The method ``AssignProperties`` takes as input the name of the element block
to which properties are being assigned, the element type for elements in the
block, the material model name, and any element fabrication properties. For
``ElasticLink2D2`` elements, the area ``A`` is the only fabrication property.

The next step is to specify the boundary conditions :math:`u_{1x}=u_{1y}=0`, and :math:`u_{12y}=0`:

.. code:: python

   step = V.StaticStep()
   step.PrescribedBC(1, (X,Y))
   step.PrescribedBC(12, Y)
   step.PrescribedBC(ALL, Z)

The point forces at nodes 3, 5, 7, 9, 11 are defined by:

.. code:: python

   step.ConcentratedLoad((3,5,9,11), Y, -10)
   step.ConcentratedLoad(7, Y, -16)

Finally, the unknown displacements are determined by solving the model and the model results are written to an ExodusII output file

   step.run()
   V.WriteResults()

The undeformed and deformed plots, generated by `ParaView <http://www.paraview.org>`__, are shown below

.. figure:: TrussExample1.png
   :align: center

   Undeformed and deformed plots of the truss. The deformed plots show contours
   of :math:`y` displacement and magnitude of the reaction forces.
