#!/usr/bin/env python
import sys
from pyfem2 import *
import matplotlib.pyplot as plt
El = PlaneStrainQuad4BBar
El = PlaneStrainQuad4Reduced

mu = 1.
Nu = .499
E = 2. * mu * (1. + Nu)

def LinearSolution(ax=None, solver=None):
    V = FiniteElementModel(jobid='VolumeLocking.Linear')
    V.GenesisMesh('QuarterCylinderQuad4.g')
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=E, Nu=Nu)
    V.InitialTemperature(ALL, 90)
    V.Temperature(ALL, 30)
    V.AssignProperties('ElementBlock1', El, 'Material-1', t=1)
    V.PrescribedBC('Nodeset-200', X)
    V.PrescribedBC('Nodeset-201', Y)
    # Pressure on inside face
    V.Pressure('Surface-1', 1.)
    V.Solve(solver=solver)
    V.WriteResults()
    ax = V.Plot2D(deformed=1, color='b' if solver is None else 'r',
                  linestyle='-.', ax=ax,
                  label='Linear + {0}'.format(solver),
                  xlim=(-.2,5), ylim=(-.2, 5))
    if solver is None:
        return ax
        ax = V.Plot2D(color='g',
                  linestyle='-', ax=ax,
                  label='undeformed',
                  xlim=(-.2,5), ylim=(-.2, 5))
    return ax

def WriteAnalyticSolution(ax=None):
    mesh = Mesh(filename='QuarterCylinderQuad4.g')
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
#ax = WriteAnalyticSolution(ax)
#ax = LinearSolution(ax)
ax = LinearSolution(ax, solver=NEWTON)
plt.legend()
plt.show()
