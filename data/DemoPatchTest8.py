from pyfem2 import *
nodtab = [[1, 0.0, 0.0],
          [2, 0.12, 0.0],
          [3, 0.24, 0.0],
          [4, 0.24, 0.06],
          [5, 0.24, 0.12],
          [6, 0.12, 0.12],
          [7, 0.0, 0.12],
          [8, 0.0, 0.06],
          [9, 0.04, 0.02],
          [10, 0.11, 0.025],
          [11, 0.18, 0.03],
          [12, 0.17, 0.055],
          [13, 0.16, 0.08],
          [14, 0.12, 0.08],
          [15, 0.08, 0.08],
          [16, 0.06, 0.05],
          [17, 0.02, 0.01],
          [18, 0.21, 0.015],
          [19, 0.2, 0.1],
          [20, 0.04, 0.1],]
eletab = [[1, 1, 3, 11, 9, 2, 18, 10, 17],
          [2, 3, 5, 13, 11, 4, 19, 12, 18],
          [3, 5, 7, 15, 13, 6, 20, 14, 19],
          [4, 7, 1, 9, 15, 8, 17, 16, 20],
          [5, 9, 11, 13, 15, 10, 12, 14, 16],]

def RunModel():
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    V = FiniteElementModel(mesh=mesh)
    mat_1 = Material('Material-1', elastic={'E': 1000000.0, 'Nu': 0.25})
    V.ElementBlock('ElementBlock-1', (1, 2, 3, 4, 5))
    V.AssignProperties('ElementBlock-1', PlaneStressQuad8, mat_1)
    step = V.StaticStep()
    step.PrescribedBC(1, X, 0.0)
    step.PrescribedBC(1, Y, 0.0)
    step.PrescribedBC(2, X, 0.00012)
    step.PrescribedBC(2, Y, 6e-05)
    step.PrescribedBC(3, X, 0.00024)
    step.PrescribedBC(3, Y, 0.00012)
    step.PrescribedBC(4, X, 0.00027)
    step.PrescribedBC(4, Y, 0.00018)
    step.PrescribedBC(5, X, 0.0003)
    step.PrescribedBC(5, Y, 0.00024)
    step.PrescribedBC(6, X, 0.00018)
    step.PrescribedBC(6, Y, 0.00018)
    step.PrescribedBC(7, X, 6e-05)
    step.PrescribedBC(7, Y, 0.00012)
    step.PrescribedBC(8, X, 3e-05)
    step.PrescribedBC(8, Y, 6e-05)
    step.run()

if __name__ == '__main__':
    RunModel()