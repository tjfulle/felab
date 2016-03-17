from conf import *
from pyfem2 import *

def test_quad4_plane_strain():
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'data/EC4SFP1.inp'))
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1e6, Nu=.25)
    V.AssignProperties('EALL', PlaneStrainQuad4, 'Material-1', t=.001)
    V.PrescribedBC(10, (X,Y), 0.)
    V.PrescribedBC(20, X, .24e-3)
    V.PrescribedBC(20, Y, .12e-3)
    V.PrescribedBC(30, X,  .3e-3)
    V.PrescribedBC(30, Y, .24e-3)
    V.PrescribedBC(40, X, .06e-3)
    V.PrescribedBC(40, Y, .12e-3)
    V.Solve()
    # Average stress must be 1600 in x and y
    step = V.steps.last
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
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'data/EC4SFP1.inp'))
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1e6, Nu=.25)
    V.AssignProperties('EALL', PlaneStressQuad4, 'Material-1', t=.001)
    V.PrescribedBC(10, (X,Y), 0.)
    V.PrescribedBC(20, X, .24e-3)
    V.PrescribedBC(20, Y, .12e-3)
    V.PrescribedBC(30, X,  .3e-3)
    V.PrescribedBC(30, Y, .24e-3)
    V.PrescribedBC(40, X, .06e-3)
    V.PrescribedBC(40, Y, .12e-3)
    V.Solve()
    step = V.steps.last
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
    V = Plane2DModel(jobid='PlaneStressTria3Patch')
    V.AbaqusMesh(filename=join(D, 'data/EC3SFP1.inp'))
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1e6, Nu=.25)
    V.AssignProperties('EALL', PlaneStressTria3, 'Material-1', t=.001)
    V.PrescribedBC(10, (X,Y), 0.)
    V.PrescribedBC(20, X, .24e-3)
    V.PrescribedBC(20, Y, .12e-3)
    V.PrescribedBC(30, X,  .3e-3)
    V.PrescribedBC(30, Y, .24e-3)
    V.PrescribedBC(40, X, .06e-3)
    V.PrescribedBC(40, Y, .12e-3)
    V.Solve()
    V.WriteResults()
    step = V.steps.last
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
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'data/EC3SFP1.inp'))
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1e6, Nu=.25)
    V.AssignProperties('EALL', PlaneStrainTria3, 'Material-1', t=.001)
    V.PrescribedBC(10, (X,Y), 0.)
    V.PrescribedBC(20, X, .24e-3)
    V.PrescribedBC(20, Y, .12e-3)
    V.PrescribedBC(30, X,  .3e-3)
    V.PrescribedBC(30, Y, .24e-3)
    V.PrescribedBC(40, X, .06e-3)
    V.PrescribedBC(40, Y, .12e-3)
    V.Solve()
    # Average stress must be 1600 in x and y
    step = V.steps.last
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
