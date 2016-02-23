.. _impOverview:

Finite element analysis stages
==============================

Finite element programs exercise three major analysis stages: preprocessing,
processing, and postprocessing. This is true of small special-purpose
programs, such as the ones developed with ``pyfem2``, and large commercial FE
programs. The general logic-flow from one analysis stage to the next is shown below

.. figure:: ProgramFlow.png
   :align: center

   Logic-flow of a linear/static FE program

Description of analysis stages
------------------------------

Preprocessing
~~~~~~~~~~~~~

Define the numerical model by assigning a mesh, materials, and boundary
conditions to the physical model problem. In this program, the
preprocessing is performed by a driver script that directly sets the
problem data

Processing
~~~~~~~~~~

Assemble and solve the system of equations describing the model problem.
Processing involves three steps:

-  Assembly of the master stiffness matrix and node force array, with
   subordinate element procedures.

-  Modification of the master stiffness matrix and node force array to
   account for prescribed boundary condtions.

-  Solution of the modified equations for unknown nodal displacements.

Postprocessing
~~~~~~~~~~~~~~

Upon executing the processing steps, nodal displacements are known. The
following postprocessing steps follow

-  Recovery of forces, including reactions
   :math:`\left(\boldsymbol{F}=\boldsymbol{K}\boldsymbol{u}\right)`.

-  Computation of internal (axial) forces and stresses in the truss
   members.

-  Plotting deflected shapes and member stress levels.
