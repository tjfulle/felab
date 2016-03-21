#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

mesh = GenesisMesh('PlateWithHoleQuad4.g')
mat = Material('Material-1')
mat.Elastic(E=100, Nu=.2)

V = FiniteElementModel(mesh=mesh, jobid='PlateWithHoleQuad4')

V.AssignProperties('ElementBlock1', PlaneStrainQuad4, mat, t=1)
V.PrescribedBC('LeftHandSide', X)
V.FixNodes('PinNode')

step = V.StaticStep()
step.SurfaceLoad(IHI, [1,0])
step.run()

V.WriteResults()

F = File('PlateWithHoleQuad4.exo')
max_p = [0., None]
max_u = [0., None]
for step in F.steps.values():
    for frame in step.frames:
        u = frame.field_outputs['U']
        for value in u.values:
            u1 = value.magnitude
            if max_u[0] < u1:
                max_u = [u1, value]

        s = frame.field_outputs['S']
        for value in s.values:
            s1 = value.max_principal
            if max_p[0] < s1:
                max_p = [s1, value]

# External and internal element numbers
xel = max_p[1].label
x = F.get_elem_coord(xel)
