.. _Material:

Material models
===============

Overview
--------

A material definition specifies the behavior of a material and supplies all
the relevant property data.

Defining the material behavior
------------------------------

Material properties are assigned the ``FiniteElementModel`` object by the ``Material`` method as:

.. code:: python

   V.Material(name)

where the string ``name`` is a name unique to the material. This command
creates an entry in the finite element space ``V``\'s ``materials`` entry for
material ``name``. Further assignments of material properties to ``name`` are
done directly on ``V.materials[name]``

Isotropic elastic properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Isotropic elastic properties are assigned to the material object by the ``Elastic`` method.  ``Elastic`` takes as input any two independent elastic constants.  For example, for a material having a Young's modulus :math:`E=1\times10^6` and Poisson's ratio :math:`\nu=.29`, the elastic properties are assigned

.. code:: python

   V.materials[name].Elastic(E=10e6, Nu=.29)

The following elastic constants are recognized by ``Elastic``:

+-----------------------+-----------+
| Property              | Key       |
+-----------------------+-----------+
| Bulk Modulus          | ``K``     |
+-----------------------+-----------+
| Shear Modulus         | ``G``     |
+-----------------------+-----------+
| Young's modulus       | ``E``     |
+-----------------------+-----------+
| Poisson's ratio       | ``Nu``    |
+-----------------------+-----------+
| First Lame parameter  | ``Lame``  |
+-----------------------+-----------+
| Second Lame parameter | ``Mu``    |
+-----------------------+-----------+

Thermal conductivity
~~~~~~~~~~~~~~~~~~~~

For thermally conductive materials, the coefficient of isotropic thermal conductivity :math:`k` is assigned by the ``IsotropicThermalConductivity`` method as:

.. code:: python

   V.materials[name].IsotropicThermalConductivity(k)

Density
~~~~~~~

The mass density :math:`\rho` of any material is assigned as

.. code:: python

   V.materials[name].Density(rho)
