.. _Loading:

Loading conditions
==================

.. _ConcentratedLoad:

Concentrated loads
------------------

.. code:: python

   step.ConcentratedLoad(region, dof[, amplitude=0.])

Distributed loads
-----------------

.. _GravityLoad:

Gravity loads
~~~~~~~~~~~~~

.. code:: python

   step.GravityLoad(region, components)

.. _DistributedLoad:

Distributed load
~~~~~~~~~~~~~~~~

.. code:: python

   step.DistributedLoad(region, components)


Surface loads
-------------

.. _SurfaceLoad:

Surface load vector
~~~~~~~~~~~~~~~~~~~

.. code:: python

   step.SurfaceLoad(surface, components)

.. _SurfaceLoadN:

Normal surface load vector
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   step.SurfaceLoadN(surface, amplitude)

.. _Pressure:

Pressure
~~~~~~~~

.. code:: python

   step.Pressure(surface, amplitude)

Heat transfer loadings
----------------------

.. _SurfaceFlux:

Surface flux
~~~~~~~~~~~~

.. code:: python

   step.SurfaceFlux(surface, qn)

.. _SurfaceConvection:

Surface convection
~~~~~~~~~~~~~~~~~~

.. code:: python

   step.SurfaceConvection(surface, Too, h)

.. _HeatGeneration:

Heat generation
~~~~~~~~~~~~~~~

.. code:: python

   step.HeatGeneration(region, amplitude)
