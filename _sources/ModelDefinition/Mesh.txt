.. _MeshDefinition:

The finite element mesh
=======================

Overview
--------

The finite element mesh

- is a collection of nodes and elements defining the geometry of the problem;
- is the fundamental object on which the finite element model is built; and
- is represented as an instance of the ``Mesh`` class in ``pyfem2``.

.. _NodeDefinition:

Nodes
-----

A node is a coordinate location in space where degrees of freedom are
defined. A node definition consists of:

-  a node label, and

-  the node coordinates.

Defining nodes
~~~~~~~~~~~~~~

A node label is a positive integer by which the node is identified. Node
labels are also referred to as external node numbers. Node labels do not
need to be numbered continuous. Individual node labels are consistently
identified by symbols or . All external communication uses external node
labels.

The node coordinates are defined as their position in a regular
Cartesian space. The coordinate system is shown in Figure
[fig:nodandel.coords].

.. figure:: Node1.png
   :align: center

   Node coordinates in Cartesian space.

Nodes are defined in a two-dimensional list where is a sublist defining
the node label and nodal coordinates of the n\ :math:`^{\rm th}` node.
For example, the node definition for node at coordinate location
:math:`P_n` in Figure [fig:nodandel.coords] is

.. code:: python

    node_n = [n, xn, yn, zn]

The collection of nodes is

.. code:: python

    nodtab = [node_1, node_2, ..., node_n, ..., node_N]

.. _ElementDefinition:

Elements
--------

An element is a mathematical relation that defines how the degrees of
freedom of a node relate to the next. An element definition consists of

-  an element type,

-  an element label, and

-  list of node labels forming the element.

The element type ID is a unique ID defined in ``pyfem``. See :ref:`ElementLibrary` for description of each element type.

Defining elements
~~~~~~~~~~~~~~~~~

An element label is a positive integer by which the element is
identified. Element labels do not need to be numbered continuously.
Individual element labels are consistently identified by symbols or .
All external communication uses external element labels.

The node labels forming the element refer to the external node labels.
Node labels must be ordered as required by the element type.

Elements are defined in a two-dimensional list where is a sublist
defining the element label, element type, and external node labels of
nodes forming element . For example, the definition for element in
Figure [fig:nodandel.elem1] is

.. figure:: Element1.png
   :align: center

   Element e

.. code:: python

    elem_e = [e, eletyp_e, i, j]

and the collection of elements is

.. code:: python

    eletab = [elem_1, elem_2, ..., elem_e, ..., elem_M]

Creating a mesh
---------------

The ``Mesh`` object can be created in several ways, as will be illustrated with the simple single element mesh of a unit square centered at the origin:

.. image:: ../apidoc/unitsquare1.png
   :scale: 35
   :align: center

Method 1
~~~~~~~~

Create the mesh by specifying ``nodtab`` and ``eletab``:

.. code:: python

   nodtab = [[10, -.5, -.5], [20, .5, -.5], [30, .5, .5], [40, -.5, .5]]
   eletab = [[100, 10, 20, 30, 40]]
   V.Mesh(nodtab=nodtab, eletab=eletab)

Method 2
~~~~~~~~

Defining a mesh manually with the ``nodtab`` and ``eletab`` is impractical for all but the simplest problems.  ``pyfem2`` allows reading in meshes created in third-party preprocessors.  To read an `ExodusII <http://prod.sandia.gov/techlib/access-control.cgi/1992/922137.pdf>`__ formatted mesh:

.. code:: python

   V.Mesh(filename='unitsquare.g')

``Mesh`` recognizes ExodusII and `vtk <www.vtk.org>`__ mesh files.  The convenience methods ``V.GenesisMesh`` and ``V.VTKMesh`` is equivalent to ``V.Mesh(filename=...)``.

Method 3
~~~~~~~~

Directly specify the nodal coordinates ``p`` and element connectivity array ``t``:

.. code:: python

   p = array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
   t = array([[0, 1, 2, 3]])
   V.Mesh(p=p, t=t)

If ``p`` and ``t`` are specified, node labels are assigned to nodes in the
order they appear in the nodal coordinates array ``p`` and take values of
``0`` to ``n-1``, where ``n`` is the length of ``p``. The element
connectivity array ``t`` must use this numbering scheme. This method is
useful for creating meshes using triangulations created by the
`pydistmesh <https://github.com/bfroehle/pydistmesh>`__ module.

Method 4
~~~~~~~~

Regular rectilinear meshes with linear quadratic elements can be created by

.. code:: python

    nx, ny = 100, 10
    lx, ly = 10, 1
    V.Mesh.RectlinearMesh((nx, ny), (lx, ly))

where ``nx`` and ``ny`` are the numer of nodes in the :math:`x` and :math:`y` coordinate directions and ``lx`` and ``ly`` are the respective lengths.

Mutual exclusivity of creation methods 1-3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``nodtab/eletab``, ``filename``, and ``p/t`` arguments are mutually exclusive.
