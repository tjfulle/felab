.. _Material:

Material models
===============

Overview
--------

A material definition specifies the behavior of a material and supplies all
the relevant property data.

Defining the material behavior
------------------------------

Material properties and responses are defined the ``Material`` object:

.. code:: python

   mat = Material(name)

where the string ``name`` is a name unique to the material. This command
creates a material ``name``. Further assignments of material properties to
``mat`` are done directly on ``mat``.

Isotropic elastic properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Isotropic elastic properties are assigned to the material object by the ``Elastic`` method.  ``Elastic`` takes as input any two independent elastic constants.  For example, for a material having a Young's modulus :math:`E=1\times10^6` and Poisson's ratio :math:`\nu=.29`, the elastic properties are assigned as

.. code:: python

   mat.Elastic(E=10e6, Nu=.29)

Alternatively, the elastic material could have also been created directly as:

.. code::

   mat = Material(name, elastic={'E':10e6, 'Nu':.29})

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

Hyperelastic properties
~~~~~~~~~~~~~~~~~~~~~~~

Hyperelastic Neo Hookean properties are assigned to the material object by the ``NeoHooke`` method.  ``NeoHooke`` takes as input any two independent elastic constants.  For example, for a material having a Young's modulus :math:`E=1\times10^6` and Poisson's ratio :math:`\nu=.29`, the elastic properties are assigned as

.. code:: python

   mat.NeoHooke(E=10e6, Nu=.29)

Alternatively, the elastic material could have also been created directly as:

.. code::

   mat = Material(name, neo_hooke={'E':10e6, 'Nu':.29})

Thermal conductivity
~~~~~~~~~~~~~~~~~~~~

For thermally conductive materials, the coefficient of isotropic thermal conductivity :math:`k` is assigned by the ``IsotropicThermalConductivity`` method as:

.. code:: python

   mat.IsotropicThermalConductivity(k)

or

.. code:: python

   mat = Material(name, thermal_conductivity=k)

Density
~~~~~~~

The mass density :math:`\rho` of any material is assigned as

.. code:: python

   mat.Density(rho)

or

.. code:: python

   mat = Material(name, density=rho)
