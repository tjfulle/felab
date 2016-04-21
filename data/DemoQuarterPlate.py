#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

V = FiniteElementModel(jobid='PlateWithHoleQuad4QuarterSym')
V.GenesisMesh('PlateWithHoleQuad4QuarterSym.g')
mat = V.Material('Material-1')
mat.Elastic(E=100, Nu=.2)
V.AssignProperties('', PlaneStrainQuad4, mat, t=1)
V.PrescribedBC('SymYZ', X)
V.PrescribedBC('SymXZ', Y)
V.InitialTemperature(ALL, 60)

step = V.StaticStep('Step-1')
step.SurfaceLoad('RightHandSide', [1,0])
step.run()

V.WriteResults()

F = File('PlateWithHoleQuad4QuarterSym.exo')
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
if not os.getenv('NOGRAPHICS'):
    #print(max_p[0])
    V.Plot2D(deformed=1)
