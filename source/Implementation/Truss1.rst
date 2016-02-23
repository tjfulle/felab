.. _TrussModel:

TrussModel: A Complete Truss Program
====================================

References
----------

-  `IFEM Chapter 2, “The Direct Stiffness Method I” <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch02.d/IFEM.Ch02.pdf>`__

-  `IFEM Chapter 3, “The Direct Stiffness Method II” <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch03.d/IFEM.Ch03.pdf>`__

Overview
--------

``TrussModel``:

-  models truss members as two-node elastic bars having constant
   cross-sectional area and Young’s modulus,

-  supports truss members defined in 1, 2, or 3 dimensions,

-  uses the direct stiffness method to assemble and solve the global
   system of equations,

-  supports prescribed displacements and nodal forces,

-  uses linear solver to solve the global system of equations, and

- writes results to `ExodusII <prod.sandia.gov/techlib/access-control.cgi/1992/922137.pdf>`__ formatted files.

.. figure:: IntroTruss.png
   :align: center

   Finite element analysis of a truss-like structure

Program limitations
~~~~~~~~~~~~~~~~~~~

-  valid only for truss-like structures whose members are two-node
   elastic bars,

-  does not support distributed loads,

-  does not support multi-freedom constraints, and

-  sacrifices speed and efficiency for clarity.

Element types
-------------

supports the following element types:

- :ref:`Link1D2` : one-dimensional two-node elastic bar.

- :ref:`Link2D2` : two-dimensional two-node elastic bar.

- :ref:`Link3D2` : three-dimensional two-node elastic bar.

The TrussProblem program
------------------------

``TrussModel`` is a finite element program that inherits from the
:ref:`apiFiniteElementModel` base class.  As such, it need only define the ``init`` and ``Solve`` methods.

The entire ``TrussModel`` script is found in ``pyfem2/fem2.py`` and is repeated here:

.. code:: python

   class TrussModel(FiniteElementModel):

       def init(self):
           # Request allocation of field variables
           self.request_output_variable('U', VECTOR, NODE)
           self.request_output_variable('R', VECTOR, NODE)
           self.request_output_variable('S', SCALAR, ELEMENT)
           self.request_output_variable('P', SCALAR, ELEMENT)

       def Solve(self):
           """Solves the finite element problem """
           # active DOF set dynamically
           self.active_dof = range(self.elements[0].ndof)
           self.validate(LinknD2, one=True)
           # Assemble the global stiffness and force
           K = self.assemble_global_stiffness()
           F, Q = self.assemble_global_force()
           Kbc, Fbc = self.apply_bc(K, F+Q)
           try:
               self.dofs = solve(Kbc, Fbc)
           except LinAlgError:
               raise RuntimeError('attempting to solve under constrained system')

           # Total force, including reaction, and reaction
           Ft = dot(K, self.dofs)
           R = Ft - F - Q

           # reshape R to be the same shape as coord
           u = self.dofs.reshape(self.numnod, -1)
           R = R.reshape(self.numnod, -1)
           p = self.internal_forces(u)
           s = self.stresses(p)
           self.snapshot(U=u, R=R, P=p, S=s)

       def internal_forces(self, u):
           p = zeros(self.numele)
           for el in self.elements:
               iel = self.mesh.elemap[el.label]
               p[iel] = el.internal_force(u[el.nodes])
           return p

       def stresses(self, p):
           return p / array([el.A for el in self.elements])

Postprocessing
~~~~~~~~~~~~~~

With nodal determined by the ``Solve`` method, postprocessing can begin. In
this program, postprocessing consists of determining the element internal
forces and stresses. The method ``TrussModel.internal_forces`` method computes
the axial internal forces of truss members.  ``TrussModel.internal_forces`` relies on the subordinate method ``internal_force`` of the ``Link`` type elements to compute the internal force (described in :ref:`apiElementLib1`).

Writing FE results
~~~~~~~~~~~~~~~~~~

The method ``WriteResults`` writes the results to a ExodusII finite element database. ExodusII files and are viewable in several commercial and open source
visualization products, including the open source
`ParaView <http://www.paraview.org>`__.  ``WriteResults`` is invoked as

.. code:: python

   V.WriteFEResults(filename)

where ``filename`` is the name of the file.  The extension '.exo' is appended to ``filename`` if it is not already present.


Example model drivers
---------------------

Example 1
~~~~~~~~~

Consider the space truss shown

.. figure:: SpaceTruss1.jpeg
   :align: center

   Example space truss

The model driver for this problem is

.. code:: python

   from pyfem2 import *
   V = TrussModel()
   nodtab = [[1, 72, 0, 0], [2, 0, 36, 0], [3, 0, 36, 72], [4, 0, 0, -48]]
   eletab = [[1, 1, 2], [2, 1, 3], [3, 1, 4]]
   V.Mesh(nodtab=nodtab, eletab=eletab)
   V.Material('Material-1')
   V.materials['Material-1'].Elastic(E=10e6, Nu=.29)
   V.ElementBlock('ElementBlock1', ALL)
   A = [.302, .729, .187]
   V.AssignProperties('ElementBlock1', Link3D2, 'Material-1', A=A)

   # Boundary conditions
   V.PrescribedBC(1, Y, 0)
   V.PrescribedBC((2,3,4), (X,Y,Z), 0)

   # Concentrated force in 'z' direction on node 1
   V.ConcentratedLoad(1, Z, -1000)

   V.Solve()

Example 2
~~~~~~~~~

Consider the bridge shown

.. figure:: Bridge.png
   :align: center

   Six-bay bridge plane truss. (a) physical problem; (b) finite element idealization. (Credit: `Ch. 21 <http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch21.d/IFEM.Ch21.pdf>`__ of Prof. Felippa’s introductory finite element course materials)

The user model definition script is

.. code:: python

    nodtab = [[1,0,0,0], [2,10,5,0], [3,10,0,0], [4,20,8,0], [5,20,0,0],
              [6,30,9,0], [7,30,0,0], [8,40,8,0], [9,40,0,0], [10,50,5,0],
              [11,50,0,0],[12,60,0,0]]
    eletab = [[1,1,3], [2,3,5], [3,5,7], [4,7,9], [5,9,11], [6,11,12],
              [7,1,2], [8,2,4], [9,4,6], [10,6,8], [11,8,10], [12,10,12],
              [13,2,3], [14,4,5], [15,6,7], [16,8,9], [17,10,11], [18,2,5],
              [19,4,7], [20,7,8], [21,9,10]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1000, Nu=.29)
    V.ElementBlock('ElementBlock1', ALL)
    Abot, Atop, Abat, Adia = 2, 10, 3, 1
    A = [Abot, Abot, Abot, Abot, Abot, Abot,
         Atop, Atop, Atop, Atop, Atop, Atop,
         Abat, Abat, Abat, Abat, Abat,
         Adia, Adia, Adia, Adia]
    V.AssignProperties('ElementBlock1', Link3D2, 'Material-1', A=A)
    V.ConcentratedLoad((3,5,9,11), Y, -10)
    V.ConcentratedLoad(7, Y, -16)
    V.PinNodes(1)
    V.PrescribedBC(ALL, Z, 0)
    V.PrescribedBC(12, Y, 0)
    V.Solve()

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

The undeformed and deformed plots as generated by ParaView:

.. figure:: TrussExample1.png
   :align: center

   Undeformed and deformed plots of the truss. The deformed plots show
   contours of :math:`y` displacement and magnitude of the reaction forces.

Example 3
~~~~~~~~~

.. figure:: TrussExample2.jpg
   :align: center

   Truss example 3

.. code:: python

    E, A = 70e9, 5 * .01 * .01
    nodtab = [[1, 0, 0], [2, 3, 4], [3, 0, 4]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=70e9, Nu=.29)
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', Link2D2, 'Material-1', A=A)
    V.PrescribedBC(1, X, -.05)
    V.PrescribedBC((2,3), (X,Y), 0)
    V.ConcentratedLoad(1, Y, 1000e3)
    V.Solve()

The nodal displacements are

.. code:: python

    [[-0.05     0.08828]
     [ 0.       0.     ]
     [ 0.       0.     ]]
