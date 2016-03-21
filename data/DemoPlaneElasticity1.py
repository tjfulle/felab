#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

V = FiniteElementModel(jobid='Plane1')
V.GenesisMesh('PlateWithHoleQuad4.g')

mat = V.Material('Material-1')
mat.Elastic(E=10e6, Nu=.29)

V.AssignProperties('ElementBlock1', PlaneStrainQuad4, mat, t=1)

step = V.StaticStep()
step.PrescribedBC('LeftHandSide', X, 0.)
step.PrescribedBC('PinNode', Y, 0.)
step.PrescribedBC('RightHandSide', X, .1)

step.run()

V.WriteResults()

if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1, colorby='Ux')
