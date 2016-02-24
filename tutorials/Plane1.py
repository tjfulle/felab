#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

V = Plane2DModel()
V.GenesisMesh('../meshes/PlateWithHoleQuad4.g')

V.Material('Material-1')
V.materials['Material-1'].Elastic(E=10e6, Nu=.29)

V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)

V.PrescribedBC('LeftHandSide', X, 0.)
V.PrescribedBC('PinNode', Y, 0.)
V.PrescribedBC('RightHandSide', X, .1)

V.Solve()

V.WriteResults('Plane1.exo')

V.Plot2D(show=1, deformed=1, colorby='Ux')
