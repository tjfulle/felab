
.. _SelectiveReducedIntegration:

Selective reduced integration elements
======================================

.. raw:: html

	 <h2> References </h2>

Overview
--------

Selective-reduced integration elements:

- form the element stiffness by integrating the deviatoric and isotropic parts of the element stiffness separately;
- reduce volumetric locking in linear and quadratic elements; and
- are the "fully integrated" elements in many production codes (see, for example, the `Abaqus Theory Guide, Ch. 3.2.4 <http://abaqus.software.polimi.it/v6.14/books/stm/default.htm>`__).


Isotropic-deviatoric split
--------------------------

The stress tensor :math:`\sigma_{ij}` (actually, any second-order tensor) can
be split additively in to isotropic and deviatoric parts:

.. math::

   \sigma_{ij} = \sigma_{ij}^{\rm iso} + \sigma_{ij}^{\rm dev}

where

.. math::

   \sigma_{ij}^{\rm iso} = \frac{\sigma_{ii}}{3}\delta_{ij} = -p\delta_{ij}, \quad
   \delta_{ij} = \begin{cases}1 & i=j\\0 & i\ne j\end{cases}

:math:`p` is the the pressure, and

.. math::

   \sigma_{ij}^{\rm dev} = \sigma_{ij} - \sigma_{ij}^{\rm dev}

Summation conventions have been used wherein repeated subscripts imply summation from 1-3.

Substituting this form of the stress tensor in to the linear elastic
constitutive equation and the finite element interpolation functions into the
weak form of the momentum equation, it can be shown that the element stiffness
matrix is

.. math::

   \left[k^e\right] = \int_{\Omega^e}
       \left[B\right]^T\left(\left[D\right] - \left[\overline{D}\right]\right)\left[B\right] dV +
   \int_{\Omega^e}
       \left[B\right]^T\left[\overline{D}\right]\left[B\right] dV

where :math:`\left[D\right]` is the material elastic stiffness and
:math:`\left[\overline{D}\right]` is

.. math::

   \overline{D}_{ij} = \frac{1}{n} \delta_i \delta_{k} D_{kj}

where :math:`n` is the number of element dimensions and

.. math::

   \delta_i = \begin{cases}
       [1, 1, 0, 0]^T & n = 2 \\
       [1, 1, 1, 0, 0, 0]^T & n = 3 \\
   \end{cases}

Selective reduced integration
-----------------------------

In a selectively reduced integration element, the first volume integral
appearing in the expression of the element stiffness is evaluated using full
integration and the second using reduced integration.
