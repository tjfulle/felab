.. _ReducedIntegration:

Reduced integration elements
============================

.. raw:: html

   <h2> References </h2>

- `Applied Mechanics of Solids, Chapter 8 <http://solidmechanics.org/Text/Chapter8_6/Chapter8_6.php#Sect8_6_2>`__

Overview
--------

Reduced integration elements:

- form the element stiffness by integrating at fewer integration points than the standard;
- reduce the computational footprint of an element;
- reduce volumetric locking in quadratic elements; and
- suffer from "hourglassing" in linear elements.

Reduced integration points
--------------------------

Reduced integration elements enjoy widespread adoption in most production
finite element codes. The concept behind reduced integration elements is
straight forward: form the element stiffness using an integration scheme that
is one order less accurate than standard. Typical choices for the number of
reduced integration points for various element types are

+-------------------------+------------------+---------------------------------------+
|Element type             | Number of nodes  | Number of reduced integration points  |
+-------------------------+------------------+---------------------------------------+
| Linear triangle         | 3                | 1                                     |
+-------------------------+------------------+---------------------------------------+
| Quadratic triangle      | 6                | 3                                     |
+-------------------------+------------------+---------------------------------------+
| Bilinear quadrilateral  | 4                | 1                                     |
+-------------------------+------------------+---------------------------------------+
+ Quadratic quadrilateral | 8                | 4                                     |
+-------------------------+------------------+---------------------------------------+

Reduced integration elements and volumetric locking
---------------------------------------------------

For the quadratic quadralateral element, reduced integration resolves
volumetric locking and even improves solution accuracy.

Hourglassing
------------

Linear reduced integration elements suffer from what is known as
"hourglassing" (for obvious reasons), as shown:

.. figure:: VolumeLocking_HG.png
   :align: center

   Hourglassing in bilinear quadrateral elements. The deformation is scaled by
   .001 in this figure (otherwise the hourglass deformation dominates).

Hourglassing is caused by weakly constrained deformation modes that lead to a
nealy singular stiffness matrix. Hourglassing must be controlled by an
:ref:`hour glass control <HourglassControl>` scheme.
