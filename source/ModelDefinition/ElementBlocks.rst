.. _ElementBlocks:

Element blocks
==============

Overview
--------

An element block is a group elements of the same type (basic geometry and
number of nodes).  All elements in the finite element model must be assigned to one and only one element block.

Defining element blocks
-----------------------

Creation of an element block requires the specification of a name and elements.

Use one of the following options to define an element block named ``name``:

- define an element block explicitly from a list of elements

.. code:: python

   elements = (1, 2, 3, ...)
   V.ElementBlock(name, elements)

the element numbers are the external element labels.

- assign all elements to an element block

.. code:: python

   V.ElementBlock(name, ALL)
