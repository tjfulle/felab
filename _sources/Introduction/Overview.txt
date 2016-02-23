Overview
========

``pyfem2`` is a Python [1]_ package providing modules, classes, and functions
for setting up and solving finite element models. It is implemented using
object oriented programming techniques. The intended audience is students who
have a familiarity with the finite element method and are interested in learning
its computational implementation.

This document serves a user’s guide to the ``pyfem2`` package and an
introduction to the computational implementation of the finite element method.  The recommended starting point is to work through the tutorials in :ref:`Tutorials` to gain a familiarity with the common ``pyfem2`` objects.  Then the reader can move on the :ref:`Implementation`.

.. [1]
   Python is an open source programming environment that comes standard
   on most operating systems and is becoming ubiqutous in computational
   sciences. However, the version of Python packaged on most operating
   systems does not come with many of the components needed for
   scientific computing (modules such as
   `numpy <http://www.numpy.org>`__ and
   `scipy <http://scipy.org/index.html>`__). Building and installing all
   of the necessary components is tedious and error prone. I recommend
   downloading and installing the
   `Anaconda <https://www.continuum.io/downloads>`__ Python distribution
   from `Continuum Analytics <https://continuum.io>`__. All source
   examples were written using Anaconda Python3.5.
