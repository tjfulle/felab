import matplotlib.pyplot as plt
from pyfem2 import *

El = PlaneStrainQuad4Reduced

E = 100
Nu = 0
V = FiniteElementModel('Test')
V.Mesh(p=[[-1,-1],[1,-1],[1,1],[-1,1]], t=[[0,1,2,3]])
V.ElementBlock('All', ALL)
V.Material('Material-1')
V.materials['Material-1'].Elastic(E=E, Nu=Nu)
V.AssignProperties('All', El, 'Material-1', t=1, hourglass_control=True)

V.PrescribedBC((0,3), X)
V.PrescribedBC((0,), Y)

step = V.StaticStep()#solver=NEWTON)
step.ConcentratedLoad(1, X, 5.)
step.ConcentratedLoad(2, X, 5.)
step.run()

V.Plot2D(deformed=True)
#print(V.dofs.reshape(-1,2))
#print(V.steps.last.frames[-1].field_outputs['S'].data[0,:,0])
plt.show()
