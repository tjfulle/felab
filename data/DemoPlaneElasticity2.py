#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

# Create the model
V = Plane2DModel()

# Read mesh from file
V.GenesisMesh('PlateWithHoleTria3.g')

V.Material('Material-1')
V.materials['Material-1'].Elastic(E=10e6, Nu=.29)

V.AssignProperties('ElementBlock1', PlaneStrainTria3, 'Material-1', t=1)

V.PrescribedBC('LeftHandSide', X, 0.)
V.PrescribedBC('BottomLeft', Y, 0)
V.PrescribedBC('RightHandSide', X, .1)

# Solve for the unknown degrees of freedom
V.Solve()

V.WriteResults('Plane2.exo')

if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1)

#PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())
