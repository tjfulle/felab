#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

V = Plane2DModel(jobid='QuadBeam')
V.RectilinearMesh((10, 2), (10, 2))
V.Material('Material-1')
V.materials['Material-1'].Elastic(E=20000, Nu=0.)
V.ElementBlock('ElementBlock1', ALL)
V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)
V.FixDOF(ILO)
V.ConcentratedLoad(IHI, Y, -10)
V.Solve()
V.WriteResults()
if not os.environ.get('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1)
