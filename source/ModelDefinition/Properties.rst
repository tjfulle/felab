.. _Properties:

Element Properties
==================

Overview
--------

Elements in an element block must be assigned material and (optionally) element fabrication properties.

Assigning element properties
----------------------------

Assignment of material and element fabrication properties requires the specification of an element block name, element type, material model name, and element fabrication properties:

.. code:: python

   V.AssignProperties(block_name, elem_type, material_name, **elem_fab)

For example, to assign elements in element block ``ElementBlock-1`` to have element type ``PlaneStrainQuad4``, material ``Material-1``, and thickness ``t=1`` use:

.. code:: python

   V.AssignProperties('ElementBlock-1', PlaneStrainQuad4, 'Material-1', t=1)

Or, to assign elements in element block ``HeatElements`` to have element type ``DiffusiveHeatTransfer2D3`` and material ``Material-1`` use:

.. code:: python

   V.AssignProperties('HeatElements', DiffusiveHeatTransfer2D3, 'Material-1')

A note on fabrication properties
--------------------------------

Element fabrication properties are required by some, but not all element types.  Element fabrication properties are properties defining the construction of an element that cannot be determined by the analysis code.  For example, plane elements require the element thickness ``t`` and beam elements require the area ``A`` and moment of inertia ``Izz``.

Element fabrication properties are defined as keyword arguments to the ``AssignProperties`` method.  Element fabrication properties can be specified as a single scalar to be applied to all elements in the block or as a list of properties - one property per element.
