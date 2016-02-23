#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

# Create the model
V = HeatTransfer2DModel()

# Read mesh from file
V.GenesisMesh('../meshes/PlateWithHoleTria3Fine.g')

# Create a material and define the thermal conductivity
V.Material('Material-1')
V.materials['Material-1'].IsotropicThermalConductivity(12)

# Define an alement block of diffusive heat transfer elements with material mat
V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')

# Fix temperatures on left and right edge
V.PrescribedBC(ILO, T, 200)
V.PrescribedBC(IHI, T, 50)

# Define surface flux on bottome edge of domain
V.SurfaceFlux(JLO, 2000)

# Define surface convection on top edge of domain
Too, h = 25, 250
V.SurfaceConvection(JHI, Too, h)

# Define a function specifying the heat generation
def fun(x):
    return 1000. / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
V.HeatGeneration(ALL, fun)

# Solve for the unknown degrees of freedom
V.Solve()

V.WriteResults('Heat1.exo')

if not os.environ.get('NOGRAPHICS'):
    V.Plot2D(colorby=V.dofs.flatten(), show=1)

#PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())
