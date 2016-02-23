
The finite element model
========================

Overview
--------

A finite element model

- is an object that stores model and history and data;
- assembles and solves the linear system representing a physical system; and
- is derived from the ``FiniteElementModel`` base class.

Finite element model classes
----------------------------

N-dimensional truss
~~~~~~~~~~~~~~~~~~~

A truss model is built from two-node elastic bar elements. In two dimensions,
a system made up of elastic bar elements is referred to as a **plane truss**.
In three dimensions, the system is referred to as a **space truss**.

.. code:: python

   V = TrussModel(mesh)

Applicable elements
...................

- :ref:`Link1D2`
- :ref:`Link2D2`
- :ref:`Link3D2`

Active degrees of freedom
.........................

- ``X``, ``Y`` (two-dimensions), ``Z`` (three-dimensions)

Valid degree of freedom assignments
...................................

- :ref:`PrescribedBC`
- :ref:`FixNodes`
- :ref:`PinNodes`

Valid loading conditions
........................

- :ref:`ConcentratedLoad`

Plane beam-column truss
~~~~~~~~~~~~~~~~~~~~~~~

The two-dimensional plane beam-column truss model is built from combinations of two-node elastic bar and Euler-Bernouli beam elements.

.. code:: python

   V = PlaneBeamColumnTrussModel(mesh)

Applicable elements
...................

- :ref:`Link2D2`
- :ref:`BeamColumn2D`

Active degrees of freedom
.........................

- ``X``, ``Y``, ``TZ``

Valid degree of freedom assignments
...................................

- :ref:`FixNodes`
- :ref:`PinNodes`
- :ref:`PrescribedBC`

Valid loading conditions
........................

- :ref:`ConcentratedLoad`

Plane elasticity
~~~~~~~~~~~~~~~~

Two-dimensional plane elasticity, in the form of plane stress and plane strain, is represented by the ``Plane2DModel`` object:

.. code:: python

   V = Plane2DModel(mesh)

Applicable elements
...................

- :ref:`PlaneStrainQuad4`
- :ref:`PlaneStrainQuad4Reduced`
- :ref:`PlaneStressQuad4`
- :ref:`PlaneStrainQuad8`
- :ref:`PlaneStressQuad8`

Active degrees of freedom
.........................

- ``X``, ``Y``

Valid degree of freedom assignments
...................................

- :ref:`FixNodes`
- :ref:`PinNodes`
- :ref:`PrescribedBC`

Valid loading conditions
........................

- :ref:`ConcentratedLoad`
- :ref:`DistributedLoad`
- :ref:`SurfaceLoad`
- :ref:`SurfaceLoadN`
- :ref:`GravityLoad`
- :ref:`Pressure`

Two-dimensional diffusive heat transfer
---------------------------------------

Two-dimensional diffusive heat transfer is modeled by the ``HeatTransfer2DModel`` object:

.. code:: python

   V = HeatTransfer2DModel(mesh)

Applicable elements
~~~~~~~~~~~~~~~~~~~

- :ref:`DiffusiveHeatTransfer2D3`

Active degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~

- ``T``

Valid degree of freedom assignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- :ref:`InitialTemperature`
- :ref:`PrescribedBC`

Valid loading conditions
~~~~~~~~~~~~~~~~~~~~~~~~

- :ref:`SurfaceFlux`
- :ref:`SurfaceConvection`
- :ref:`HeatGeneration`
