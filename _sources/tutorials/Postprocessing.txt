.. _PostProcessing:

Post processing
===============

Objective
---------

This example demonstrates the use of the ``pyfem2`` scripting interface to the ``ExodusII`` output database.

Problem description
-------------------

Determine the maximum principal stress in the plate with hole and compare with the analytic solution.

.. code:: python

   from pyfem2 import *

   V = FiniteElementModel(jobid='PlateQuad')
   V.GenesisMesh('PlateWithHoleQuad4QuarterSym.g')

   mat = V.Material('Material-1', elastic={'E':100, 'Nu':.2})
   V.AssignProperties('ElementBlock1', PlaneStrainQuad4, mat, t=1)
   V.PrescribedBC('LeftHandSide', X)
   V.FixNodes('PinNode')

   # Create a load step, apply a traction to the right hand side, and run it
   step = V.StaticStep()
   step.SurfaceLoad(IHI, [1,0])
   step.run()
   V.WriteResults()

   # Open the output file
   F = File('PlateQuad.exo')

   # Loop through each load step, determining the maximum principal stress
   max_u = [0., None]
   max_p = [0., None]
   for step in F.steps:
       u = step.field_outputs['U']
       for value in u.values:
           u1 = value.magnitude
           if max_u[0] < u1:
               max_u = [u1, value]

       s = step.field_outputs['S']
       for value in s.values:
           s1 = value.max_principal
           if max_p[0] < s1:
               max_p = [s1, value]

   # External and internal element numbers
   xel = max_p[1].label
   x = F.get_elem_coord(xel)

How does it work?
-----------------

The complete code can be found in ``pyfem2/data/Postprocessing.py``. We now examine the preceding program in detail, beginning with the ``F = File(...)`` line.

The line

.. code:: python

   # Open the output file
   F = File('PlateQuad.exo')

opens and loads the previously written ExodusII file and instantiates a ``ExodusII`` file object.

The statements

.. code:: python

   max_u = [0., None]
   max_p = [0., None]

create containers to hold the maximum displacement magnitude and maximum principal stress.

.. code:: python

   for step in F.steps:
       u = step.field_outputs['U']
       for value in u.values:
           u1 = value.magnitude
           if max_u[0] < u1:
               max_u = [u1, value]

       s = step.field_outputs['S']
       for value in s.values:
           s1 = value.max_principal
           if max_p[0] < s1:
               max_p = [s1, value]

Get the element label and coordinates of the element having the maximum principal stress.

.. code:: python

   xel = max_p[1].label
   x = F.get_elem_coord(xel)
