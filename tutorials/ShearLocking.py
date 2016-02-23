#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

mu = 10000.
nu = 0.
E = 2. * mu * (1. + nu)

def WriteAnalyticSolution():
    mesh = Mesh(filename='../meshes/SlenderBeamQuad4.g')
    a = 0.15
    b = 1.
    P = 2 * a * b
    L = 10.
    d = -P * L ** 3 / (2 * E * a ** 3 * b)
    II = b * (2. * a) ** 3 / 12.
    u = zeros_like(mesh.coord)
    for (i, x) in enumerate(mesh.coord):
        x1, x2 = x
        x2 -= a
        u[i,0] = P / (2. * E * II) * x1 ** 2 * x2 + \
            nu * P / (6 * E * II) * x2 ** 3 - \
            P / (6 * II * mu) * x2 ** 3 - \
            (P * L ** 2 / (2 * E * II) - P * a ** 2 / (2 * II * mu)) * x2
        u[i,1] = -nu * P / (2 * E * II) * x1 * x2 ** 2 - \
            P / (6 * E * II) * x1 ** 3 + \
            P * L ** 2 / (2 * E * II) * x1 - P * L ** 3/(3 * E * II)
    mesh.PutNodalSolution('ShearLocking.Analytic', u)
    return

def LinearSolution():
    V = Plane2DModel()
    V.GenesisMesh('../meshes/SlenderBeamQuad4.g')
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=E, Nu=nu)
    V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)
    V.PrescribedBC('Nodeset-200', (X,Y))
    V.PrescribedBC('Nodeset-201', X)
    V.SurfaceLoad('Surface-300', [0, -1])
    V.Solve()
    V.WriteResults('ShearLocking.Linear')

WriteAnalyticSolution()
LinearSolution()
