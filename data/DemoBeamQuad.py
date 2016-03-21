#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

mesh = RectilinearMesh2D(nx=10, ny=2, lx=10, ly=2)
mat = Material('Material-1', elastic={'E':20000, 'Nu':0})

V = FiniteElementModel(mesh=mesh, jobid='QuadBeam')
V.ElementBlock('ElementBlock1', ALL)
V.AssignProperties('ElementBlock1', PlaneStrainQuad4, mat, t=1)

step = V.StaticStep()
step.FixDOF(ILO)
step = V.StaticStep()
step.ConcentratedLoad(IHI, Y, -10)
step.run()
V.WriteResults()
if not os.environ.get('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1)
