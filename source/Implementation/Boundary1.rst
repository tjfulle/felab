Boundary conditions and concentrated loads
==========================================

Overview
--------

Boundary conditions applied to the finite element model are represented internally as ``doftags`` and ``dofvals``.

Internal representation of boundary conditions and concentrated loads
---------------------------------------------------------------------

Internally, boundary condition data for every node in the model are
stored in two arrays:

- ``doftags`` specifies the degree of freedom information; and

- ``dofvals`` specifies the magnitudes of prescribed displacement/concentrated
  force boundary conditions.

The arrays are populated dynamically as boundary conditions and concentrated loads are assigned to the finite element model.

The arrays ``doftags`` and ``dofvals`` are best explained by
example. Consider the space truss:

.. image:: SpaceTruss1.jpeg
   :align: center
   :scale: 50

The boundary conditions and concentrated loads are defined by

.. code:: python

    # Boundary conditions
    V.PinNodes((2,3,4))
    V.PrescribedBC(1, Y, 0)
    # Concentrated force in 'z' direction on node 1
    V.ConcentratedLoad(1, Z, -1000)

The concentrated force of *0* applied to the ``X`` DOF of node 1 was not
explicitly specified. Internally, nodes on which displacements or point forces
are not explicitly prescribed are assumed to have a point force of 0 applied.

The corresponding ``doftags`` and ``dofvals`` arrays are:

.. code:: python

   >>> print doftags
   [[0  1  0  0  0  0  0]
    [1  1  1  0  0  0  0]
    [1  1  1  0  0  0  0]
    [1  1  1  0  0  0  0]]
   >>> print dofvals
   [[0  0  -1000  0  0  0  0]
    [0  0      0  0  0  0  0]
    [0  0      0  0  0  0  0]
    [0  0      0  0  0  0  0]]

``doftag[i,j]`` is the *j*\ :sup:`th` degree of freedom tag for the *i*\
:sup:`th` node. ``doftags[i,j]=1`` indicates that degree of freedom *j* of
node *i* is prescribed while ``doftags[i,j]=0`` indicates that the
corresponding force is known.

The magnitudes of the prescribed degrees of freedom or force are transferred
to the ``dofvals`` array. ``dofvals[i,j]`` is the magnitude corresponding to
the prescribed condition on the *j*\ :sup:`th` degree of freedom of node
*i*.

.. note::

   Every node for which a force/displacement boundary condition is not
   explicitly prescribed is assigned a nodal force of *0*.
