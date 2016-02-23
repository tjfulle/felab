.. _InternalNodesAndElements:

Nodes and elements
==================

References
----------

- :ref:`Defining nodes <NodeDefinition>`
- :ref:`Defining elements <ElementDefinition>`

Overview
--------

Nodes and elements are assigned to the finite element model in one of several ways.  Internally, node and element data are stored in several arrays and mappings.

Internal representation of nodes and elements
---------------------------------------------

Nodes
~~~~~

Internally, nodes are assigned integer IDs ranging from *0* to *n*, where *n*
is the number of nodes in the problem. Node IDs are assigned in the order they
appear in the user defined ``nodtab`` table. Node coordinates and internal
node IDs are stored separately in

- the nodal coordinate array ``coord`` with ``coord[i]`` being the coordinates of
  of the *i*\ :sup:`th` node; and

- the node ID mapping ``nodmap`` where ``inode=nodmap[xnode]`` is the internal
  node ID of the external node label ``xnode``.

In ``pyfem2``, internal node IDs are consistently identified by the symbol
``inode``.

Elements
~~~~~~~~

Like nodes, elements are assigned integer IDs ranging from *0* to *m*, where
*m* is the number of elements in the problem. Internal element IDs are
assigned in the order that they are assigned to element blocks (not the order
they appear in the ``eletab`` table). Element node lists, internal element
IDs, and element object lists are stored separately in:

- the element connectivity array ``elecon``, where ``elecon[i,j]`` is the
  internal node ID of the *j*\ :sup:`th` node of element *i*;

- the element ID mapping ``elemap``, where ``ielem=elemap[xelem]`` is the
  internal element ID of the external element label ``xelem``; and

- the element object list ``elements``, where ``elements[ielem]`` is the
  instantiated element object of element ``ielem``.

In ``pyfem2``, internal element IDs are consistently identified by the symbol
``xelem`` (or ``xe``)).

Generating the internal node and element representations
--------------------------------------------------------

The method :ref:`parse_nod_and_elem_tables <parse_nodes_and_elements>`
generates the internal node and element representations. The ``coord`` and
``nodmap`` objects are created directly from the ``nodtab`` table. Since
internal element IDs are generated in the order elements are assigned to
element blocks, they are not generated in ``parse_nod_and_elem_tables``.
Instead, a mapping from external element numbers to the corresponding internal
node numbers ``elemap1`` is created.

Example
~~~~~~~

Consider the space truss shown

.. figure:: SpaceTruss1.jpeg
   :align: center
   :scale: 60

   Example space truss

The node and element tables are

.. code:: python

    nodtab = [[1, 72,  0,   0],
              [2,  0, 36,   0],
              [3,  0, 36,  72],
              [4,  0,  0, -48]]
    eletab = [[1, 1, 2],
              [2, 1, 3],
              [3, 1, 4]]

The corresponding ``coord`` array and ``nodtab`` mapping are:

.. code:: python

   >>> print coord
   [[72.  0.   0.],
    [ 0. 36.   0.],
    [ 0. 36.  72.],
    [ 0.  0. -48.]]
   >>> print nodmap
   {1:0, 2:1, 3:2, 4:3}

As noted, the mapping from external to internal element numbers cannot be created until after elements have been assigned to element blocks.  The mapping ``elemap1`` is:

.. code:: python

   >>> print elemap1
   {1: [0, 1],
    2: [0, 2],
    3: [0, 3]}
