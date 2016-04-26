from pyfem2 import *
nodtab = [[1, 0.0, 0.0],
          [2, 0.24, 0.0],
          [23, 0.24, 0.12],
          [4, 0.0, 0.12],
          [15, 0.04, 0.02],
          [6, 0.18, 0.03],
          [7, 0.16, 0.08],
          [100, 0.08, 0.08],]
eletab = [[1, 1, 2, 6],
          [2, 1, 6, 15],
          [3, 4, 1, 15],
          [4, 4, 15, 100],
          [5, 2, 23, 7],
          [6, 2, 7, 6],
          [7, 23, 4, 100],
          [8, 23, 100, 7],
          [9, 15, 6, 7],
          [10, 15, 7, 100],]

def RunModel():
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    V = FiniteElementModel(mesh=mesh)
    mat_1 = Material('Material-1', elastic={'Nu': 0.25, 'E': 1000000.0})
    V.ElementBlock('ElementBlock-1', (1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
    V.AssignProperties('ElementBlock-1', PlaneStressTria3, mat_1)
    step = V.StaticStep()
    step.PrescribedBC(1, X, 0.0)
    step.PrescribedBC(1, Y, 0.0)
    step.PrescribedBC(2, X, 0.00024)
    step.PrescribedBC(2, Y, 0.00012)
    step.PrescribedBC(4, X, 6e-05)
    step.PrescribedBC(4, Y, 0.00012)
    step.PrescribedBC(23, X, 0.0003)
    step.PrescribedBC(23, Y, 0.00024)
    step.run()

if __name__ == '__main__':
    RunModel()