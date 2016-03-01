
.. _SelectiveReducedIntegration:

Selective reduced integration elements
======================================

.. raw:: html

	 <h2> References </h2>

- `The Finite Element Method, Ch 4.4 <https://books.google.com/books?id=cHH2n_qBK0IC&pg=PA221&lpg=PA221&dq=hughes+selective+reduced+integration+google+book&source=bl&ots=vrupu77YTR&sig=Yps9ajb6yVA_lOnqbkF2C1XwhWI&hl=en&sa=X&ved=0ahUKEwiv_NSjgaDLAhVI8x4KHSf0BaMQ6AEIUjAJ#v=onepage&q=hughes%20selective%20reduced%20integration%20google%20book&f=false>`__

Overview
--------

Selective-reduced integration elements:

- form the element stiffness by integrating the deviatoric and isotropic parts of the element stiffness separately;
- reduce volumetric locking in linear and quadratic elements; and
- are the "fully integrated" elements in many production codes (see, for example, the `Abaqus Theory Guide, Ch. 3.2.4 <http://abaqus.software.polimi.it/v6.14/books/stm/default.htm>`__).


Selective reduced integration
-----------------------------

The element stiffness :math:`\left[k^e\right]` is given by

.. math::

   \boldsymbol{k}_{ij} =
   \int_{\Omega^e}\boldsymbol{B}_i^T\boldsymbol{D}\boldsymbol{B}_j\, d\Omega

The matrix of material properties can split additively in to distortional and dilatational parts as

.. math::

   \boldsymbol{D} = \overline{\boldsymbol{D}} + \hat{\boldsymbol{D}}

where

- **Plane strain**

  .. math::

     \overline{\boldsymbol{D}} = \mu \begin{bmatrix}
        2 &   &   &   \\
          & 2 &   &   \\
          &   & 2 &   \\
          &   &   & 1
     \end{bmatrix}
     \quad
     \hat{\boldsymbol{D}} = \lambda \begin{bmatrix}
        1 & 1 & 1 & 0 \\
        1 & 1 & 1 & 0 \\
        1 & 1 & 1 & 0 \\
        0 & 0 & 0 & 0
     \end{bmatrix}

- **Three dimensions**

  .. math::

     \overline{\boldsymbol{D}} = \mu \begin{bmatrix}
        2 &   &   &   &   &  \\
          & 2 &   &   &   &  \\
          &   & 2 &   &   &  \\
          &   &   & 1 &   &  \\
          &   &   &   & 1 &  \\
          &   &   &   &   & 1
     \end{bmatrix}
     \quad
     \hat{\boldsymbol{D}} = \lambda \begin{bmatrix}
        1 & 1 & 1 & 0 & 0 & 0 \\
        1 & 1 & 1 & 0 & 0 & 0 \\
        1 & 1 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0
     \end{bmatrix}

where :math:`\mu` is the shear modulus :math:`\mu=\frac{E}{2(1+\nu)}` and
:math:`\lambda` is the first Lame parameter
:math:`\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}`

.. note::

   The additive split shown is valid for strain tensors stored as :math:`\boldsymbol{\epsilon}=\begin{bmatrix} \epsilon_{xx} & \epsilon_{yy} & \epsilon_{zz} & 2\epsilon_{xy}\end{bmatrix}^T` for plane strain and :math:`\boldsymbol{\epsilon}=\begin{bmatrix} \epsilon_{xx} & \epsilon_{yy} & \epsilon_{zz} & 2\epsilon_{xy} & 2\epsilon_{yz} & 2\epsilon_{xz}\end{bmatrix}^T` for three dimensions.

With this decomposition, the element stiffness can be expressed as

.. math::

   \boldsymbol{k}_{ij} = \overline{\boldsymbol{k}}_{ij} + \hat{\boldsymbol{k}}_{ij}

where

.. math::

   \overline{\boldsymbol{k}}_{ij} =
   \int_{\Omega^e}\boldsymbol{B}_i^T\overline{\boldsymbol{D}}\boldsymbol{B}_j\, d\Omega

.. math::

   \hat{\boldsymbol{k}}_{ij} =
   \int_{\Omega^e}\boldsymbol{B}_i^T\hat{\boldsymbol{D}}\boldsymbol{B}_j\, d\Omega

For nearly incompressible materials, since :math:`\lambda/\mu >> 1`, the
numerical values of :math:`\hat{\boldsymbol{k}}` tend to be much greater than
those of :math:`\overline{\boldsymbol{k}}`. Consequently,
:math:`\hat{\boldsymbol{k}}` dominates the element response and leads to
volumetric locking if nonzero volumetric strains are computed at any of the
element's integration points. A simple method of alleviating this tendency is
to reduce the order of numerical quadrature used to evaluate
:math:`\hat{\boldsymbol{k}}`. This technique is known as **selective reduced
integration**.
