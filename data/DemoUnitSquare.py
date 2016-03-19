import os
from pyfem2 import *

V = FiniteElementModel()
V.UnitSquareMesh((1,1))
V.ElementBlock('All', ALL)
V.Material('Mat', elastic={'E': 1000, 'Nu': 0})
V.AssignProperties('All', PlaneStressQuad4, 'Mat')
V.FixNodes(ILO)

step = V.StaticStep()
step.PrescribedBC(IHI, X, .1)
step.run()
if not os.getenv('NOGRAPHICS'):
    V.Plot2D(deformed=1, show=1)

step = V.StaticStep()
step.PrescribedBC(IHI, X, 0)
step.run()

if not os.getenv('NOGRAPHICS'):
    V.Plot2D(deformed=1, show=1)