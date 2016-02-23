.. _Loading:

Loading conditions
==================

.. _ConcentratedLoad:

Concentrated loads
------------------

.. code:: python

   V.ConcentratedLoad(region, dof[, amplitude=0.])

Distributed loads
-----------------

.. _GravityLoad:

Gravity loads
~~~~~~~~~~~~~

.. code:: python

   V.GravityLoad(region, components)

.. _DistributedLoad:

Distributed load
~~~~~~~~~~~~~~~~

.. code:: python

   V.DistributedLoad(region, components)


Surface loads
-------------

.. _SurfaceLoad:

Surface load vector
~~~~~~~~~~~~~~~~~~~

.. code:: python

   V.SurfaceLoad(surface, components)

.. _SurfaceLoadN:

Normal surface load vector
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   V.SurfaceLoadN(surface, amplitude)

.. _Pressure:

Pressure
~~~~~~~~

.. code:: python

   V.Pressure(surface, amplitude)

Heat transfer loadings
----------------------

.. _SurfaceFlux:

Surface flux
~~~~~~~~~~~~

.. code:: python

   V.SurfaceFlux(surface, qn)

.. _SurfaceConvection:

Surface convection
~~~~~~~~~~~~~~~~~~

.. code:: python

   V.SurfaceConvection(surface, Too, h)

.. _HeatGeneration:

Heat generation
~~~~~~~~~~~~~~~

.. code:: python

   V.HeatGeneration(region, amplitude)
