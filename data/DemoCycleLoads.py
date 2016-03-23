from pyfem2 import *

mesh = UnitSquareMesh(nx=1, ny=1)
mat = Material('Mat-1', elastic={'E': 500, 'Nu': 0})

V = FiniteElementModel(mesh=mesh)
V.ElementBlock('Block-1', ALL)
V.AssignProperties('Block-1', PlaneStressQuad4, mat)

V.FixNodes(ILO)

step = V.StaticStep(solver=NEWTON)
step.PrescribedBC(IHI, X, .1)
step.run()
if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=True)

step = V.StaticStep(solver=NEWTON)
step.PrescribedBC(IHI, X, 0)
step.run()
if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=True)

step = V.StaticStep(solver=NEWTON)
step.RemoveBC(IHI, X)
step.SurfaceLoad(IHI, [100, 0])
step.run()
if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=True)
