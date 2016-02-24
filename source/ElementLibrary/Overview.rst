
Overview
========

``pyfem2`` provides several elements in one- and two-dimensions.

Heat transfer elements
----------------------

.. _DiffusiveHeatTransfer2D3:

``DiffusiveHeatTransfer2D3``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``DiffusiveHeatTransfer2D3`` element is a three-node plane element that is appropriate for diffusive heat transfer problems in 2D.

- **Active degrees of freedom:** ``T``
- **Valid loads:** :ref:`HeatGeneration`, :ref:`SurfaceFlux`, :ref:`SurfaceConvection`

Truss elements
--------------

The truss element theory is described in :ref:`BarElements`.

.. _Link1D2:

``ElasticLink1D2``
~~~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X``
- **Valid loads:** :ref:`ConcentratedLoad`

.. _Link2D2:

``ElasticLink2D2``
~~~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y``
- **Valid loads:** :ref:`ConcentratedLoad`

.. _Link3D2:

``ElasticLink3D2``
~~~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y,Z``
- **Valid loads:** :ref:`ConcentratedLoad`

Beam-column elements
--------------------

.. _BeamColumn2D:

``BeamColumn2D``
~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y,TZ``
- **Valid loads:** :ref:`ConcentratedLoad`

Stress-displacement elements
----------------------------

.. _PlaneStrainQuad4:

PlaneStrainQuad4
~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y``
- **Valid loads:** :ref:`ConcentratedLoad`, :ref:`GravityLoad`, :ref:`DistributedLoad`, :ref:`SurfaceLoad`, :ref:`SurfaceLoadN`, :ref:`Pressure`

.. _PlaneStrainQuad4Reduced:

PlaneStrainQuad4Reduced
~~~~~~~~~~~~~~~~~~~~~~~

Uses a reduced integration formulation.

- **Active degrees of freedom:** ``X,Y``
- **Valid loads:** :ref:`ConcentratedLoad`, :ref:`GravityLoad`, :ref:`DistributedLoad`, :ref:`SurfaceLoad`, :ref:`SurfaceLoadN`, :ref:`Pressure`

.. _PlaneStressQuad4:

PlaneStressQuad4
~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y``
- **Valid loads:** :ref:`ConcentratedLoad`, :ref:`GravityLoad`, :ref:`DistributedLoad`, :ref:`SurfaceLoad`, :ref:`SurfaceLoadN`, :ref:`Pressure`

.. _PlaneStrainQuad8:

PlaneStrainQuad8
~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y``
- **Valid loads:** :ref:`ConcentratedLoad`, :ref:`GravityLoad`, :ref:`DistributedLoad`, :ref:`SurfaceLoad`, :ref:`SurfaceLoadN`, :ref:`Pressure`

.. _PlaneStressQuad8:

PlaneStressQuad8
~~~~~~~~~~~~~~~~

- **Active degrees of freedom:** ``X,Y``
- **Valid loads:** :ref:`ConcentratedLoad`, :ref:`GravityLoad`, :ref:`DistributedLoad`, :ref:`SurfaceLoad`, :ref:`SurfaceLoadN`, :ref:`Pressure`
