.. _Install:

Installation
============

Obtaining ``pyfem2``
--------------------

``pyfem2`` is a pure Python package requiring only , in addition to the
standard library. may be obtained in several ways, but the easiest is to
download and install the `Anaconda <http://www.continuum.io>`__ Python
Distribution. is also easily installed on linux using the distribution’s
package manager.

``pyfem2`` is installed by first

-  downloading the package from `Git Hub <https://github.com/tjfulle/pyfem2>`__, or

-  cloning the package from `<https://github.com/tjfulle/pyfem2>`__ with ``git``.

Installing ``pyfem2``
---------------------

After obtaining ``pyfem2``, your Python distribution must be able to find
it. There are several ways to accomplish this. Perhaps the easiest is to
move the package to a location on your file system and add that location
to your ``PYTHONPATH`` environment variable. For example, if ``pyfem2`` is downloaded
and placed in the ``Documents`` folder in your home directory:

- **Mac OS X or Linx:** add ``export PYTHONPATH=$HOME/Documents/pyfem2`` to your ``.bashrc`` (if using a variant of ``bash``) or ``setenv PYTHONPATH $HOME/Documents/pyfem2`` to your ``.cshrc`` (if using a variant of ``csh``)

-  **Windows:** create a ``PYTHONPATH`` variable that points to ``%HOMEPATH%\My Documents\pyfem2`` environment variable (``My Computer > Properties > Advanced System Settings > Environment Variables``)

Testing the installation
------------------------

Test the installation by navigating to ``$HOME/Documents/pyfem2/examples`` and executing

.. code:: shell

   python plane1.py

at a command prompt.

If the script executes without error, the installation was successful.

If ``pyfem2`` was copied to a location other than ``$HOME/Documents/pyfem2``, adjust the paths in these instructions accordingly.
