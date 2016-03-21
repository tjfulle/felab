
The finite element model
========================

Overview
--------

A finite element model

- is an instance of the ``FiniteElementModel`` class;
- stores model and history and data; and
- assembles and solves the linear system representing a physical system.

The finite element model object
-------------------------------

Finite element models are represented by the ``FiniteElementModel`` class.  ``FiniteElementModel`` is a model neutral container that sets up and stores data needed to run a finite element model.  Load steps are assigned to the model and run.


Static stress-displacement step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   step = V.StaticStep()

Applicable elements
...................

- :ref:`Link1D2`
- :ref:`Link2D2`
- :ref:`Link3D2`

Active degrees of freedom
.........................

- ``X``
- ``Y``
- ``Z``
- ``TX``
- ``TY``
- ``TZ``

Valid degree of freedom assignments
...................................

- :ref:`PrescribedBC`
- :ref:`FixNodes`
- :ref:`PinNodes`

Valid loading conditions
........................

- :ref:`ConcentratedLoad`
- :ref:`DistributedLoad`
- :ref:`SurfaceLoad`
- :ref:`SurfaceLoadN`
- :ref:`GravityLoad`
- :ref:`Pressure`

Heat transfer step
------------------

.. code:: python

   V.HeatTransferStep()

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
