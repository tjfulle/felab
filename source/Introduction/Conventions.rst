Conventions
===========

The following conventions that are used throughout ``pyfem2`` are
discussed:

-  Degrees of freedom
-  Symbolic constants

Degrees of freedom
------------------

The “degrees of freedom,” (DOF) of a physical system describe its
spatial configuration. For example, the spatial configuration of a
mechanical system could be described by the position, rotation, and
temperature of every point in a system. If the number of degrees of
freedom is finite, the model is said to be discrete, and continuous
otherwise. Because FEM is a discretization method, the number of DOF of
a FEM model is necessarily finite.

The degrees of freedom at a point are referred to as follows:

#. x-displacement

#. y-displacement

#. z-displacement

#. Rotation about the x-axis, in radians

#. Rotation about the y-axis, in radians

#. Rotation about the z-axis, in radians

#. Temperature

The x-, y-, and z-directions coincide with the global X-, Y-, and
Z-directions, respectively. Not every degree of freedom is activated for
every point in space. For example, two-dimensional bodies subject only
to displacements in the x- and y-directions use only the degrees of
freedom 1 and 2.

Symbolic constants
------------------

Symbolic constants are defined in ``pyfem2`` to identify common objects, e.g.,
degree of freedom labels (among others). Symbolic constants are written in all
capital letters. Symbolic constants must be imported in to a script before
using them with the following statement:

.. code:: python

    from pyfem2.constants import *

Importing the entire ``pyfem2`` namespace with

.. code:: python

    from pyfem2 import *

also imports the symbolic constants.

Symbolic constants representing degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``X``, ``Y``, ``Z`` represent the :math:`x`, :math:`y`, and :math:`z` Cartesian coordinate directions, respectively.

- ``TX``, ``TY``, ``TZ`` represent the rotation :math:`theta` about the :math:`x`, :math:`y`, and :math:`z` Cartesian axes, respectively.

- ``T`` represents the temperature.

Symbolic constants representing geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``ALL`` represents all geometric entities, e.g., all nodes, all elements, etc.

- ``ILO``, ``IHI``, ``JLO``, ``JHI`` correspond to the :math:`x` (``I``) and a :math:`y` (``J``) coordinate directions and the identifiers ``LO`` and ``HI`` to the corresponding low and high boundaries.

- ``BOUNDARY`` correspond to entities on the finite element domain boundary.

- ``S1``, ``S2``, ... ``S10`` represent the :math:`{\rm i}^{\rm{th}}` side of an element.

Symbolic constants representing data types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``SCALAR`` represents scalar data.

- ``VECTOR`` represents vector data.

- ``SYMTENSOR`` represents symmetric tensor data.

Symbolic constants representing field positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``NODE`` represents the nodal field position.

- ``ELEMENT`` represents the element field position.

- ``ELEMENT_CENTROID`` represents the element centroid field position.

- ``INTEGRATION_POINT`` represents the integration point field position.
