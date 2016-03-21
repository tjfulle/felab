#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

# Create the model
V = FiniteElementModel(jobid='Heat1')

# Read mesh from file
V.GenesisMesh('PlateWithHoleTria3Fine.g')

# Create a material and define the thermal conductivity
mat = V.Material('Material-1')
mat.IsotropicThermalConductivity(12)

# Define an alement block of diffusive heat transfer elements with material mat
V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, mat)

# Fix temperatures on left and right edge
step = V.HeatTransferStep()
step.PrescribedBC(ILO, T, 200)
step.PrescribedBC(IHI, T, 50)

# Define surface flux on bottome edge of domain
step.SurfaceFlux(JLO, 2000)

# Define surface convection on top edge of domain
Too, h = 25, 250
step.SurfaceConvection(JHI, Too, h)

# Define a function specifying the heat generation
def fun(x):
    return 1000. / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
step.HeatGeneration(ALL, fun)

# Solve for the unknown degrees of freedom
step.run()

V.WriteResults()

if not os.environ.get('NOGRAPHICS'):
    V.Plot2D(colorby=step.dofs.flatten(), show=1)

#PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())
