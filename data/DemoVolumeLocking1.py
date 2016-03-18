#!/usr/bin/env python
import os
import sys
from pyfem2 import *
import matplotlib.pyplot as plt

def Runall(ax=None):
    mu = 1.
    plt.clf()
    plt.cla()
    for (i, Nu) in enumerate((0., .2, .45, .499)):
        E = 2. * mu * (1. + Nu)

        filename = '../../pyfem2-docs/_images/VolumeLocking_%d.png'%i

        # Analytic solution
        mesh = Mesh(filename='QuarterCylinderQuad4.g')
        p = 1.
        a, b = mesh.coord[0, 1], mesh.coord[-1, 0]
        u = zeros_like(mesh.coord)
        for (i, x) in enumerate(mesh.coord):
            r = sqrt(x[0] ** 2 + x[1] ** 2)
            term1 = (1. + Nu) * a ** 2 * b ** 2 * p / (E * (b ** 2 - a ** 2))
            term2 = 1. / r + (1. - 2. * Nu) * r / b ** 2
            ur = term1 * term2
            u[i, :] = ur * x[:] / r
        mesh.PutNodalSolution('VolumeLocking.Analytic', u)
        ax = mesh.Plot2D(xy=u+mesh.coord, color='orange',
                         label=r'Analytic, $\nu=%g$'%Nu, weight=8)

        # Linear finite element solution
        V = FiniteElementModel()
        V.GenesisMesh('QuarterCylinderQuad4.g')
        V.Material('Material-1')
        V.materials['Material-1'].Elastic(E=E, Nu=Nu)
        V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1')
        V.PrescribedBC('Nodeset-200', X)
        V.PrescribedBC('Nodeset-201', Y)

        step = V.StaticStep()
        # Pressure on inside face
        step.Pressure('Surface-1', 1.)
        step.run()
        V.Plot2D(deformed=1, color='b', linestyle='-.',
                 label=r'Linear, $\nu=%g$'%Nu, ax=ax, filename=filename,
                 xlim=(-.2,5), ylim=(-.2, 5))

Runall()
