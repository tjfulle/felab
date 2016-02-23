Defining a model
================

Overview
--------

A finite element model in ``pyfem2`` is created by:

- creating a ``FiniteElementModel`` object;
- assigning a finite element mesh to the model object;
- assigning elements to element blocks;
- assigning boundary conditions; and
- [optionally] applying some form of external loading.

The ``FiniteElementModel`` object
---------------------------------

``FiniteElementModel`` is the base class from which all other finite element models are derived.  The ``FiniteElementModel`` object builds the finite element model on a finite element mesh.  Specific subclasses of ``FiniteElementModel`` are described in :ref:`Modeling`.

The finite element mesh
-----------------------

The mesh of a model is defined by elements and their nodes. The rules and methods for defining nodes and elements are described in :ref:`NodeDefinition`, :ref:`ElementDefinition`, and :ref:`MeshDefinition`.

Element blocks
--------------

An element block is a group of elements having the same element type and
material definition. Every element must be assigned to one, and only one,
element block.  Creation of element blocks is described in :ref:`ElementBlocks`.

Boundary conditions
-------------------

Zero and/or non-zero value boundary conditions must be imposed on one or more model degrees of freedom.  Specification of boundary conditions is described in :ref:`Boundary`.

External loading
----------------

Usually, some form of external loading is applied to the model, e.g., concentrated or distributed loads.  Assignment of loading conditions is described in :ref:`Loading`.
