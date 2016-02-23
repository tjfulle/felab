.. _Sets:

Entity sets
===========

Overview
--------

It is convenient to group nodes, elements, and surfaces in model so that
boundary conditions and loads can be applied to entire groups, rather than
individual entitities.

.. _NodeSets:

Node sets
---------

A node set is a group of nodes to which boundary conditions and concentrated loads can be applied.  A node set is created by:

.. code:: python

   V.NodeSet(name, nset)

where ``name`` is a unique character string by which the set is referenced and ``nset`` is one of the following:

- a single external node label:

  .. code:: python

     V.NodeSet(name, 1)

- a list of external node labels:

  .. code:: python

     V.NodeSet(name, (1,2,3,...))

- a symbolic constant:

  .. code:: python

     V.NodeSet(name, CONSTANT)

  where ``CONSTANT`` is one of the following symbolic constants:

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

Element sets
------------

An element set is a group of elements to which distributed loads can be applied.  A element set is created by:

.. code:: python

   V.ElementSet(name, elset)

where ``name`` is a unique character string by which the set is referenced and ``elset`` is one of the following:

- a single external element label:

  .. code:: python

     V.ElementSet(name, 1)

- a list of external element labels:

  .. code:: python

     V.ElementSet(name, (1,2,3,...))

- a symbolic constant:

  .. code:: python

     V.ElementSet(name, CONSTANT)

  where ``CONSTANT`` is one of the following symbolic constants:

  +---------------+-------------------------------------------------------+
  | Constant      | Description                                           |
  +---------------+-------------------------------------------------------+
  | ``ALL``       | All of the elements in the mesh                       |
  +---------------+-------------------------------------------------------+

Surfaces
--------

A surface is a group of element/face pairs to which surface loads can be applied.  A surface is created by:

.. code:: python

   V.Surface(name, surface)

where ``name`` is a unique character string by which the surface is referenced and ``surface`` is one of the following:

- a single external element label/face ID pair:

  .. code:: python

     V.Surface(name, (1, S1))

  which specifies that the surface is defined by the face ``S1`` of element 1. In general, the symbolic constant ``Sn`` represents the ``n``\th face of the element.

- a list of external element labels/face IDs:

  .. code:: python

     V.Surface(name, ((1, S1), (2, S1), (3, S1), ...))

- one of the following symbolic constants:

  +---------------+-------------------------------------------------------+
  | Constant      | Description                                           |
  +---------------+-------------------------------------------------------+
  | ``BOUNDARY``  | All of the surfaces on the boundary of the mesh       |
  +---------------+-------------------------------------------------------+
  | ``ILO``       | Low boundary in the :math:`x` coordinate direction    |
  +---------------+-------------------------------------------------------+
  | ``IHI``       | High boundary in the :math:`x` coordinate direction   |
  +---------------+-------------------------------------------------------+
  | ``JLO``       | Low boundary in the :math:`y` coordinate direction    |
  +---------------+-------------------------------------------------------+
  | ``JHI``       | High boundary in the :math:`y` coordinate direction   |
  +---------------+-------------------------------------------------------+
