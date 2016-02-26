from conf import *
from pyfem2 import *

def test_quad4_plane_strain():
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'meshes/EC4SFP1.inp'))
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
    step = V.steps[-1]
    field = step.field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1600.)
        assert allclose(data[:,1], 1600.)
        assert allclose(data[:,2], 800.)
        assert allclose(data[:,3], 400.)
    field = step.field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,3], 1e-3)

def test_quad4_plane_stress():
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'meshes/EC4SFP1.inp'))
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
    step = V.steps[-1]
    field = step.field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.33333333)
        assert allclose(data[:,1], 1333.33333333)
        assert allclose(data[:,2], 400.)
    field = step.field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,2], 1e-3)

def test_tria3_plane_stress():
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'meshes/EC3SFP1.inp'))
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
    #V.Plot2D(show=1)
    V.Solve()
    V.WriteResults('PlaneStressTria3Patch.exo')
    step = V.steps[-1]
    field = step.field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.333333333), 'Wrong Sxx'
        assert allclose(data[:,1], 1333.333333333), 'Wrong Syy'
        assert allclose(data[:,2],  400.), 'Wrong Sxy'
    field = step.field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,2], 1e-3)

def test_tria3_plane_strain():
    V = Plane2DModel()
    V.AbaqusMesh(filename=join(D, 'meshes/EC3SFP1.inp'))
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
    step = V.steps[-1]
    field = step.field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1600.)
        assert allclose(data[:,1], 1600.)
        assert allclose(data[:,2], 800.)
        assert allclose(data[:,3], 400.)
    field = step.field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,3], 1e-3)
