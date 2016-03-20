from conf import *
from pyfem2 import *

def test_quad4_plane_strain():
    V = FiniteElementModel()
    V.AbaqusMesh(filename=join(D, 'data/EC4SFP1.inp'))
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.AssignProperties('EALL', PlaneStrainQuad4, mat, t=.001)
    step = V.StaticStep()
    step.PrescribedBC(10, (X,Y), 0.)
    step.PrescribedBC(20, X, .24e-3)
    step.PrescribedBC(20, Y, .12e-3)
    step.PrescribedBC(30, X,  .3e-3)
    step.PrescribedBC(30, Y, .24e-3)
    step.PrescribedBC(40, X, .06e-3)
    step.PrescribedBC(40, Y, .12e-3)
    step.run()

    # Average stress must be 1600 in x and y
    field = step.frames[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1600.)
        assert allclose(data[:,1], 1600.)
        assert allclose(data[:,2], 800.)
        assert allclose(data[:,3], 400.)
    field = step.frames[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,3], 1e-3)

def test_quad4_plane_stress():
    V = FiniteElementModel()
    V.AbaqusMesh(filename=join(D, 'data/EC4SFP1.inp'))
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.AssignProperties('EALL', PlaneStressQuad4, mat, t=.001)
    step = V.StaticStep()
    step.PrescribedBC(10, (X,Y), 0.)
    step.PrescribedBC(20, X, .24e-3)
    step.PrescribedBC(20, Y, .12e-3)
    step.PrescribedBC(30, X,  .3e-3)
    step.PrescribedBC(30, Y, .24e-3)
    step.PrescribedBC(40, X, .06e-3)
    step.PrescribedBC(40, Y, .12e-3)
    step.run()
    field = step.frames[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.33333333)
        assert allclose(data[:,1], 1333.33333333)
        assert allclose(data[:,2], 400.)
    field = step.frames[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,2], 1e-3)

def test_tria3_plane_stress():
    V = FiniteElementModel(jobid='PlaneStressTria3Patch')
    V.AbaqusMesh(filename=join(D, 'data/EC3SFP1.inp'))
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.AssignProperties('EALL', PlaneStressTria3, mat, t=.001)
    step = V.StaticStep()
    step.PrescribedBC(10, (X,Y), 0.)
    step.PrescribedBC(20, X, .24e-3)
    step.PrescribedBC(20, Y, .12e-3)
    step.PrescribedBC(30, X,  .3e-3)
    step.PrescribedBC(30, Y, .24e-3)
    step.PrescribedBC(40, X, .06e-3)
    step.PrescribedBC(40, Y, .12e-3)
    step.run()
    V.WriteResults()
    field = step.frames[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.333333333), 'Wrong Sxx'
        assert allclose(data[:,1], 1333.333333333), 'Wrong Syy'
        assert allclose(data[:,2],  400.), 'Wrong Sxy'
    field = step.frames[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,2], 1e-3)

def test_tria3_plane_strain():
    V = FiniteElementModel()
    V.AbaqusMesh(filename=join(D, 'data/EC3SFP1.inp'))
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.AssignProperties('EALL', PlaneStrainTria3, mat, t=.001)
    step = V.StaticStep()
    step.PrescribedBC(10, (X,Y), 0.)
    step.PrescribedBC(20, X, .24e-3)
    step.PrescribedBC(20, Y, .12e-3)
    step.PrescribedBC(30, X,  .3e-3)
    step.PrescribedBC(30, Y, .24e-3)
    step.PrescribedBC(40, X, .06e-3)
    step.PrescribedBC(40, Y, .12e-3)
    step.run()
    # Average stress must be 1600 in x and y
    field = step.frames[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1600.)
        assert allclose(data[:,1], 1600.)
        assert allclose(data[:,2], 800.)
        assert allclose(data[:,3], 400.)
    field = step.frames[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,3], 1e-3)
