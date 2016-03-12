#!/usr/bin/env python
from pyfem2 import *
import matplotlib.pyplot as plt

mu = 10000.
nu = 0.
E = 2. * mu * (1. + nu)

def PlaneStressBeam(ratio):
    V = Plane2DModel()
    length = 10.
    a = ratio * length
    P = 2.22 / length ** 3 * E * a ** 3
    q = P / a
    V.RectilinearMesh((10, 3), (length, 2*a))
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=E, Nu=nu)
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', PlaneStressQuad4Incompat, 'Material-1', t=1)
    V.PrescribedBC(IHI, (X,Y))
    V.SurfaceLoad(ILO, [0, -q])
    V.Solve()
    if not os.getenv('NOGRAPHICS'):
        ax = V.Plot2D(deformed=1, color='b')
    u = AnalyticSolution(V.mesh.coord, q)
    xy = V.mesh.coord + u
    if not os.getenv('NOGRAPHICS'):
        V.mesh.Plot2D(xy=xy, color='r', ax=ax, show=1)

def AnalyticSolution(coord, q):
    a = (coord[:,1].max() - coord[:,1].min()) / 2.
    L = coord[:,0].max() - coord[:,0].min()
    b = 1.
    P = 2 * q * a * b
    d = -P * L ** 3 / (2 * E * a ** 3 * b)
    II = b * (2. * a) ** 3 / 12.
    u = zeros_like(coord)
    for (i, x) in enumerate(coord):
        x1, x2 = x
        x2 -= a
        u[i,0] = P / (2. * E * II) * x1 ** 2 * x2 + \
            nu * P / (6 * E * II) * x2 ** 3 - \
            P / (6 * II * mu) * x2 ** 3 - \
            (P * L ** 2 / (2 * E * II) - P * a ** 2 / (2 * II * mu)) * x2
        u[i,1] = -nu * P / (2 * E * II) * x1 * x2 ** 2 - \
            P / (6 * E * II) * x1 ** 3 + \
            P * L ** 2 / (2 * E * II) * x1 - P * L ** 3/(3 * E * II)
    return u

# ---- Change the aspect ratio to see shear locking
ratio = .0015
PlaneStressBeam(ratio)
