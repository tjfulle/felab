
.. _PlaneElasticity1:

Plane Elasticity, 1
===================

Objective
---------

This example demonstrates the use of the ``Plane2DModel`` to solve a plane strain problem.

Problem Description
-------------------

A beam is approximated by plane strain elements with a thickness of 1 and the
material is linear elastic with constant Young's modulus :math:`E` and
Poisson's ratio :math:`\nu`. The beam is fixed at its left end and concentrated forces applied to the right.

.. image:: QuadBeam0.png
   :align: center
   :scale: 80

Model Script
------------

The problem is defined and solved in ``pyfem2`` as follows:

.. code:: python

   """
   pyfem2 tutorial demo program: Plane strain elasticity

   """

   from pyfem2 import *

   # Create the finite element model and rectilinear mesh
   V = Plane2DModel()
   V.RectilinearMesh((10, 2), (10, 2))

   # Assign element to element block
   V.ElementBlock('ElementBlock1', ALL)

   # Create a material and define the elastic properties
   V.Material('Material-1')
   V.materials['Material-1'].Elastic(E=20000, Nu=0.)
   V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)

   # Fix the left end apply loads to right
   V.FixDOF(ILO)
   V.ConcentratedLoad(IHI, Y, -10)

   # Solve the model and write results to the ExodusII database
   V.Solve()
   V.WriteResults('QuadBeam.exo')

How does it work?
-----------------

We now examine the preceding program in detail.

The first lines of the program,

.. code:: python

   from pyfem2 import *

import objects from the ``pyfem2`` namespace in to the program. The key import
 from the ``pyfem2`` library is ``Plane2DModel``.

The statement

.. code:: python

   V = Plane2DModel()

creates the finite element model.  A rectilinear mesh with 10 elements in the :math:`x` direction and 2 elements in the :math:`y` with lengths :math:`L_x=10` and :math:`L_y=10` is created by

.. code:: python

   V.RectilinearMesh((10, 2), (10, 2))

.. image:: QuadBeam1.png
   :align: center
   :scale: 80

The finite element model ``V`` requires that every element in the mesh be
assigned an element type and material constitutive relations. The assignment
occurs by grouping elements of the same type together in to element blocks and
then assigning to those element blocks material and fabrication properties.
For this problem, every element in the mesh is a ``PlaneStrainQuad4`` four-node
plane strain element (a quadralateral element with two degrees of freedom per
node). The element block, named ``ElementBlock1``, containing all of the
elements in the mesh is created by:

.. code:: python

   V.ElementBlock('ElementBlock1', ALL)

The symbol ``ALL`` informs the ``ElementBlock`` method to assign all elements
in the mesh to the block ``ElementBlock1``. Material properties are defined by
the ``Material`` method and elements in a block are assigned material and
fabrication properties by the ``AssignProperties`` method:

.. code:: python

   V.Material('Material-1')
   V.materials['Material-1'].Elastic(E=20000, Nu=0.)
   V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)

The method ``AssignProperties`` takes as input the name of the element block
to which properties are being assigned, the element type for elements in the
block, the material model name, and any element fabrication properties. For
``PlaneStrainQuad4`` elements, the thickness ``t`` is the only fabrication
property.

The next step is to specify the boundary conditions at :math:`x=0` and the concentrated nodal forces at :math:`x=10`:

.. code:: python

   V.FixDOF(ILO)
   V.ConcentratedLoad(IHI, Y, -10)

The symbols ``ILO`` and ``IHI`` correspond to the :math:`x` coordinate
direction (``I``) and the identifiers ``LO`` and ``HI`` to the corresponding
low and high boundaries.

Finally, the unknown displacements are determined by solving the model and the model results are written to an ExodusII output file

.. code:: python

   V.Solve()

Perhaps the easiest way to view results is by:

.. code:: python

   V.Plot2D(deformed=1, show=1)

With the ``deformed`` keyword, the deformed coordinates are plotted.

.. image:: QuadBeam2.png
   :align: center
   :scale: 80

The results can also be written to an ExodusII file and viewed in
`ParaView <http://www.paraview.org>`__

.. code:: python

   V.WriteResults('QuadBeam.exo')

.. image:: QuadBeam3.png
   :align: center
