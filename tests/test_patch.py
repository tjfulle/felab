from conf import *
from felab import *

def test_quad4_plane_strain():
    V = fe_model()
    V.abaqus_mesh(filename='./data/EC4SFP1.inp')
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.assign_properties('EALL', CPE4, mat, t=.001)
    stage = V.create_static_stage()
    stage.assign_prescribed_bc(10, (X,Y), 0.)
    stage.assign_prescribed_bc(20, X, .24e-3)
    stage.assign_prescribed_bc(20, Y, .12e-3)
    stage.assign_prescribed_bc(30, X,  .3e-3)
    stage.assign_prescribed_bc(30, Y, .24e-3)
    stage.assign_prescribed_bc(40, X, .06e-3)
    stage.assign_prescribed_bc(40, Y, .12e-3)
    stage.run()

    # Average stress must be 1600 in x and y
    field = stage.increments[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1600.)
        assert allclose(data[:,1], 1600.)
        assert allclose(data[:,2], 800.)
        assert allclose(data[:,3], 400.)
    field = stage.increments[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,3], 1e-3)

def test_quad4_plane_stress():
    V = fe_model()
    V.abaqus_mesh(filename='./data/EC4SFP1.inp')
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.assign_properties('EALL', CPS4, mat, t=.001)
    stage = V.create_static_stage()
    stage.assign_prescribed_bc(10, (X,Y), 0.)
    stage.assign_prescribed_bc(20, X, .24e-3)
    stage.assign_prescribed_bc(20, Y, .12e-3)
    stage.assign_prescribed_bc(30, X,  .3e-3)
    stage.assign_prescribed_bc(30, Y, .24e-3)
    stage.assign_prescribed_bc(40, X, .06e-3)
    stage.assign_prescribed_bc(40, Y, .12e-3)
    stage.run()
    field = stage.increments[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.33333333)
        assert allclose(data[:,1], 1333.33333333)
        assert allclose(data[:,2], 400.)
    field = stage.increments[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,2], 1e-3)

def test_tria3_plane_stress():
    V = fe_model(jobid='PlaneStressTria3Patch')
    V.abaqus_mesh(filename='./data/EC3SFP1.inp')
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.assign_properties('EALL', CPS3, mat, t=.001)
    stage = V.create_static_stage()
    stage.assign_prescribed_bc(10, (X,Y), 0.)
    stage.assign_prescribed_bc(20, X, .24e-3)
    stage.assign_prescribed_bc(20, Y, .12e-3)
    stage.assign_prescribed_bc(30, X,  .3e-3)
    stage.assign_prescribed_bc(30, Y, .24e-3)
    stage.assign_prescribed_bc(40, X, .06e-3)
    stage.assign_prescribed_bc(40, Y, .12e-3)
    stage.run()
    V.write_results()
    field = stage.increments[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.333333333), 'Wrong Sxx'
        assert allclose(data[:,1], 1333.333333333), 'Wrong Syy'
        assert allclose(data[:,2],  400.), 'Wrong Sxy'
    field = stage.increments[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,2], 1e-3)

def test_tria3_plane_strain():
    V = fe_model()
    V.abaqus_mesh(filename='./data/EC3SFP1.inp')
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})
    V.assign_properties('EALL', CPE3, mat, t=.001)
    stage = V.create_static_stage()
    stage.assign_prescribed_bc(10, (X,Y), 0.)
    stage.assign_prescribed_bc(20, X, .24e-3)
    stage.assign_prescribed_bc(20, Y, .12e-3)
    stage.assign_prescribed_bc(30, X,  .3e-3)
    stage.assign_prescribed_bc(30, Y, .24e-3)
    stage.assign_prescribed_bc(40, X, .06e-3)
    stage.assign_prescribed_bc(40, Y, .12e-3)
    stage.run()
    # Average stress must be 1600 in x and y
    field = stage.increments[-1].field_outputs['S']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1600.)
        assert allclose(data[:,1], 1600.)
        assert allclose(data[:,2], 800.)
        assert allclose(data[:,3], 400.)
    field = stage.increments[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3)
        assert allclose(data[:,1], 1e-3)
        assert allclose(data[:,3], 1e-3)
