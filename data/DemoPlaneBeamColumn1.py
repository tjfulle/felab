from pyfem2 import *

nodtab = [[1,-4,3], [2,0,0], [3,0,3], [4,nan,nan], [5,4,3]]
eletab = [[1,1,3], [2,3,5], [3,1,2], [4,2,3], [5,2,5]]
mesh = Mesh(nodtab=nodtab, eletab=eletab)

Ec, Em = 30000, 200000
mat1 = Material('Material-1', elastic={'E':Ec, 'Nu':.3})
mat2 = Material('Material-2', elastic={'E':Em, 'Nu':.3})

V = FiniteElementModel(mesh=mesh, jobid='PlaneBeamColumn')
V.ElementBlock('B1', (1,2))
V.ElementBlock('B2', (3,5))
V.ElementBlock('B3', (4,))
V.AssignProperties('B1', PlaneBeamColumn, mat1, A=.02, Izz=.004)
V.AssignProperties('B2', ElasticLink2D2, mat2, A=.001)
V.AssignProperties('B3', ElasticLink2D2, mat2, A=.003)

V.PrescribedBC(1, (X,Y,TZ))
V.PrescribedBC(5, Y)

step = V.StaticStep()
step.ConcentratedLoad(2, Y, 100)
step.ConcentratedLoad(5, TZ, 200)
step.ConcentratedLoad(5, X, 300)
step.run()
V.WriteResults()
