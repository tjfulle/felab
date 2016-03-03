#!/usr/bin/env python
import sys
from pyfem2 import *
import matplotlib.pyplot as plt

mu = 1.
Nu = .499
E = 2. * mu * (1. + Nu)

def LinearSolution(ax=None):
    V = Plane2DModel()
    V.GenesisMesh('../meshes/QuarterCylinderQuad4.g')
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=E, Nu=Nu)
    V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)
    V.PrescribedBC('Nodeset-200', X)
    V.PrescribedBC('Nodeset-201', Y)
    # Pressure on inside face
    V.Pressure('Surface-1', 1.)
    V.Solve()
    V.WriteResults('VolumeLocking.Linear')
    ax = V.Plot2D(deformed=1, color='b', linestyle='-.', ax=ax, label='Linear',
                  xlim=(-.2,5), ylim=(-.2, 5))
    return ax

def ReducedIntegrationSolution(ax=None):
    V = Plane2DModel()
    V.GenesisMesh('../meshes/QuarterCylinderQuad4.g')
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=E, Nu=Nu)
    V.AssignProperties('ElementBlock1', PlaneStrainQuad4Reduced, 'Material-1', t=1)
    V.PrescribedBC('Nodeset-200', X)
    V.PrescribedBC('Nodeset-201', Y)
    # Pressure on inside face
    V.Pressure('Surface-1', 1.)
    V.Solve()
    V.WriteResults('VolumeLocking.Reduced')
    ax = V.Plot2D(deformed=1, ax=ax, color='b', label='Reduced integration')
    return ax

def QuadraticSolution(ax=None):
    V = Plane2DModel()
    V.GenesisMesh('../meshes/QuarterCylinderQuad8.g')
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=E, Nu=Nu)
    V.AssignProperties('ElementBlock1', PlaneStrainQuad8, 'Material-1', t=1)
    V.PrescribedBC('Nodeset-200', X)
    V.PrescribedBC('Nodeset-201', Y)
    # Pressure on inside face
    V.Pressure('Surface-1', 1.)
    #V.SurfaceLoad("Surface-300", [0.195090322, 0.98078528])
    #V.SurfaceLoad("Surface-301", [0.555570233, 0.831469612])
    #V.SurfaceLoad("Surface-302", [0.831469612, 0.555570233])
    #V.SurfaceLoad("Surface-303", [0.98078528, 0.195090322])
    V.Solve()
    V.WriteResults('VolumeLocking.Quadratic')

def WriteAnalyticSolution(ax=None):
    mesh = Mesh(filename='../meshes/QuarterCylinderQuad4.g')
    a = mesh.coord[0, 1]
    b = mesh.coord[-1, 0]
    p = 1.
    u = zeros_like(mesh.coord)
    for (i, x) in enumerate(mesh.coord):
        r = sqrt(x[0] ** 2 + x[1] ** 2)
        term1 = (1. + Nu) * a ** 2 * b ** 2 * p / (E * (b ** 2 - a ** 2))
        term2 = 1. / r + (1. - 2. * Nu) * r / b ** 2
        ur = term1 * term2
        u[i, :] = ur * x[:] / r
    mesh.PutNodalSolution('VolumeLocking.Analytic', u)
    ax = mesh.Plot2D(xy=u+mesh.coord, ax=ax, color='orange', weight=8,
                     label='Analytic')
    return ax

ax = None
ax = WriteAnalyticSolution(ax)
ax = ReducedIntegrationSolution(ax)
ax = LinearSolution(ax)
QuadraticSolution()

if not os.environ.get('NOGRAPHICS'):
    plt.legend()
    plt.show()
