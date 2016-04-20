from conf import *
from pyfem2 import *

def test_4_node_plane_stress():
    nodtab = [[ 0, 0.0,  0.0],
              [ 1, 0.24, 0.0],
              [22, 0.24, 0.12],
              [ 3, 0.0,  0.12],
              [14, 0.04, 0.02],
              [ 5, 0.18, 0.03],
              [ 6, 0.16, 0.08],
              [99, 0.08, 0.08]]
    eletab = [[1, 0, 1, 5, 14],
              [2, 1, 22, 6, 5],
              [3, 22, 3, 99, 6],
              [4, 3, 0, 14, 99],
              [5, 14, 5, 6, 99]]
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    mat = Material('Mat-1', elastic={'E':1e6, 'Nu':.25})

    V = FiniteElementModel(mesh=mesh)
    V.ElementBlock('Block-1', ALL)
    V.AssignProperties('Block-1', PlaneStressQuad4, mat)

    # Fixed BC
    V.FixNodes(0)

    # Load step
    step = V.StaticStep()

    step.PrescribedBC(1, X, 2.4e-4)
    step.PrescribedBC(1, Y, 1.2e-4)

    step.PrescribedBC(22, X, 3.0e-4)
    step.PrescribedBC(22, Y, 2.4e-4)

    step.PrescribedBC(3, X, 6.0e-5)
    step.PrescribedBC(3, Y, 1.2e-4)

    step.run()

    s = step.frames[-1].field_outputs['S'].data
    assert allclose(s[:,:,:2], 1333.3333)
    assert allclose(s[:,:,2], 400.0)

def test_8_node_plane_stress():
    nodtab = [[ 0, 0.0,  0.0],
              [ 1, 0.12, 0.0],
              [ 2, 0.24, 0.0],
              [ 3, 0.24, 0.06],
              [ 4, 0.24, 0.12],
              [ 5, 0.12, 0.12],
              [ 6, 0.0,  0.12],
              [ 7, 0.0,  0.06],
              [ 8, 0.04, 0.02],
              [ 9, 0.11, 0.025],
              [10, 0.18, 0.03],
              [11, 0.17, 0.055],
              [12, 0.16, 0.08],
              [13, 0.12, 0.08],
              [14, 0.08, 0.08],
              [15, 0.06, 0.05],
              [16, 0.02, 0.01],
              [17, 0.21, 0.015],
              [18, 0.20, 0.10],
              [19, 0.04, 0.10]]
    eletab = [[1, 0, 2,  10, 8,  1, 17,  9, 16],
              [2, 2, 4,  12, 10, 3, 18, 11, 17],
              [3, 4, 6,  14, 12, 5, 19, 13, 18],
              [4, 6, 0,  8,  14, 7, 16, 15, 19],
              [5, 8, 10, 12, 14, 9, 11, 13, 15]]

    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    mat = Material('Mat-1', elastic={'E':1e6, 'Nu':.25})

    V = FiniteElementModel(mesh=mesh, jobid='foo')
    V.ElementBlock('Block-1', ALL)
    V.AssignProperties('Block-1', PlaneStressQuad8, mat)

    V.FixNodes(0)

    step = V.StaticStep()
    step.PrescribedBC(1, X, 1.2e-4)
    step.PrescribedBC(1, Y, 6.0e-5)
    step.PrescribedBC(2, X, 2.4e-4)
    step.PrescribedBC(2, Y, 1.2e-4)
    step.PrescribedBC(3, X, 2.7e-4)
    step.PrescribedBC(3, Y, 1.8e-4)
    step.PrescribedBC(4, X, 3.0e-4)
    step.PrescribedBC(4, Y, 2.4e-4)
    step.PrescribedBC(5, X, 1.8e-4)
    step.PrescribedBC(5, Y, 1.8e-4)
    step.PrescribedBC(6, X, 6.0e-5)
    step.PrescribedBC(6, Y, 1.2e-4)
    step.PrescribedBC(7, X, 3.0e-5)
    step.PrescribedBC(7, Y, 6.0e-5)

    step.run()

    s = step.frames[-1].field_outputs['S'].data
    assert allclose(s[:,:,:2], 1333.3333)
    assert allclose(s[:,:,2], 400.0)

def test_3_node_plane_stress():
    nodtab = [[ 0,  0.0, 0.0],
              [ 1, 0.24, 0.0],
              [22, 0.24, 0.12],
              [ 3,  0.0, 0.12],
              [14, 0.04, 0.02],
              [ 5, 0.18, 0.03],
              [ 6, 0.16, 0.08],
               [99, 0.08, 0.08]]
    eletab = [[ 1,  0,  1, 5],
              [ 2,  0,  5, 14],
              [ 3,  3,  0, 14],
              [ 4,  3, 14, 99],
              [ 5,  1, 22, 6],
              [ 6,  1,  6, 5],
              [ 7, 22,  3, 99],
              [ 8, 22, 99, 6],
              [ 9, 14,  5, 6],
              [10, 14,  6, 99]]
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    mat = Material('Mat-1', elastic={'E':1e6, 'Nu':.25})

    V = FiniteElementModel(mesh=mesh)
    V.ElementBlock('Block-1', ALL)
    V.AssignProperties('Block-1', PlaneStressTria3, mat)

    V.FixNodes(0)

    step = V.StaticStep()
    step.PrescribedBC(1,  X, 2.4e-4)
    step.PrescribedBC(1,  Y, 1.2e-4)
    step.PrescribedBC(22, X, 3.0e-4)
    step.PrescribedBC(22, Y, 2.4e-4)
    step.PrescribedBC(3,  X, 6.0e-5)
    step.PrescribedBC(3,  Y, 1.2e-4)

    step.run()

    s = step.frames[-1].field_outputs['S'].data
    assert allclose(s[:,:,:2], 1333.3333)
    assert allclose(s[:,:,2], 400.0)
