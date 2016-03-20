#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

V = FiniteElementModel(jobid='QuadBeam')
V.RectilinearMesh((10, 2), (10, 2))
mat = V.Material('Material-1')
mat.Elastic(E=20000, Nu=0.)
V.ElementBlock('ElementBlock1', ALL)
V.AssignProperties('ElementBlock1', PlaneStrainQuad4, mat.name, t=1)
V.FixDOF(ILO)

step = V.StaticStep()
step.ConcentratedLoad(IHI, Y, -10)
step.run()
V.WriteResults()
if not os.environ.get('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1)
