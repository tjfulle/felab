.. _Boundary:

Boundary conditions
===================

Overview
--------

Boundary conditions:

- are used to specify the values of basic solution variables (displacements, rotations, temperature) at nodes;

- are defined *after* elements have been assigned to element blocks.

.. _PrescribedBC:

Prescribing degrees of freedom
------------------------------

Degrees of freedom can be prescribed by:

.. code:: python

   V.PrescribedBC(region, dof[, amplitude=0])

- ``region`` is a node label, list of node labels, node set, or one of the symbolic constants:

  +---------------+-------------------------------------------------------+
  | Constant      | Description                                           |
  +---------------+-------------------------------------------------------+
  | ``ALL``       | All of the nodes in the mesh                          |
  +---------------+-------------------------------------------------------+
  | ``BOUNDARY``  | All of the nodes on the boundary of the mesh          |
  +---------------+-------------------------------------------------------+
  | ``ILO``       | Low boundary in the :math:`x` coordinate direction    |
  +---------------+-------------------------------------------------------+
  | ``IHI``       | High boundary in the :math:`x` coordinate direction   |
  +---------------+-------------------------------------------------------+
  | ``JLO``       | Low boundary in the :math:`y` coordinate direction    |
  +---------------+-------------------------------------------------------+
  | ``JHI``       | High boundary in the :math:`y` coordinate direction   |
  +---------------+-------------------------------------------------------+

- ``dof`` is the degree of freedom direction, one of the symbolic constants ``X``, ``Y``, ``Z``, ``TX``, ``TY``, ``TZ``, or ``T``.  Can be given as a list of degrees of freedom.

- the optional ``amplitude`` is the magnitude of the prescribed boundary condition.  ``amplitude`` is either a scalar constant or function of the coordinates.

For example, to fix the :math:`x` displacement of the node set ``Nodeset-1`` do:

.. code:: python

   V.PrescribedBC('Nodeset-1', X)

Or, to prescribe a temperature in the :math:`x` and :math:`y` given by :math:`2\frac{1+y}{(3+x)^2+(1+y)^2}` on the boundary of a mesh do:

.. code:: python

     def fun(x):
         return 2. * (1. + x[:,1]) / ((3. + x[:,0])**2 + (1 + x[:,1])**2)
     fun = lambda x: sqrt(x[:,0]**2 + x[:,1]**2)
     V.PrescribedBC(BOUNDARY, T, fun)

.. _PinNodes:

Pinning displacement degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Displacement degrees of freedom can be fixed by:

.. code:: python

   V.PinNodes(region)

where ``region`` is a node label, list of node labels, node set, or symbolic constant (see :ref:`PrescribedBC`.)

.. _FixNodes:

Fixing displacement degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Displacement and rotation degrees of freedom can be fixed by:

.. code:: python

   V.FixNodes(region)

where ``region`` is a node label, list of node labels, node set, or symbolic constant (see :ref:`PrescribedBC`.)

.. _InitialTemperature:

Initial temperature
-------------------

The initial temperature at nodes can be set by:

.. code:: python

   V.InitialTemperature(region, amplitude)
