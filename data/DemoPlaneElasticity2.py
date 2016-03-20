#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

# Create the model
V = FiniteElementModel(jobid='Plane2')

# Read mesh from file
V.GenesisMesh('PlateWithHoleTria3.g')

mat = V.Material('Material-1')
mat.Elastic(E=10e6, Nu=.29)

V.AssignProperties('ElementBlock1', PlaneStrainTria3, mat.name, t=1)

step = V.StaticStep()
step.PrescribedBC('LeftHandSide', X, 0.)
step.PrescribedBC('BottomLeft', Y, 0)
step.PrescribedBC('RightHandSide', X, .1)

# Solve for the unknown degrees of freedom
step.run()

V.WriteResults()

if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1)

#PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())
