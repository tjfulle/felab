import matplotlib.pyplot as plt
from pyfem2 import *

El = PlaneStrainQuad4SelectiveReduced

E = 100
Nu = 0
V = Plane2DModel('Test')
V.Mesh(p=[[-1,-1],[1,-1],[1,1],[-1,1]], t=[[0,1,2,3]])
V.ElementBlock('All', ALL)
V.Material('Material-1')
V.materials['Material-1'].Elastic(E=E, Nu=Nu)
V.AssignProperties('All', El, 'Material-1', t=1)

V.PrescribedBC((0,3), X)
V.PrescribedBC((0,), Y)
V.ConcentratedLoad(1, X, 5.)
V.ConcentratedLoad(2, X, 5.)
#V.Solve(solver=NEWTON, increments=2)
V.Solve()
V.Plot2D(deformed=True)
#print(V.dofs.reshape(-1,2))
#print(V.steps.last.frames[-1].field_outputs['S'].data[0,:,0])
plt.show()
