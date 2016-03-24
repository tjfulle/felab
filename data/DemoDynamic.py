from pyfem2 import *

mesh = UnitSquareMesh(nx=1, ny=1)
mat = Material('Mat-1', elastic={'E': 500, 'Nu': 0})
mat.Density(1.)

V = FiniteElementModel(mesh=mesh)
V.ElementBlock('Block-1', ALL)
V.AssignProperties('Block-1', PlaneStrainQuad4, mat)

step = V.DynamicStep(period=1e-6, increments=10)
step.PrescribedBC(IHI, X, .1)
step.run()
if not os.getenv('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=True)
