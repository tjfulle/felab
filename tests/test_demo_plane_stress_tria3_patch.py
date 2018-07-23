from felab import *

def PlaneStressTria3PatchDemo(plot=False):
    # READ MESH
    mesh = abaqus_mesh(filename='./data/EC3SFP1.inp')

    # CREATE MATERIAL MODEL
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})

    # CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
    V = fe_model(mesh=mesh, jobid='Tri3PlaneStressPatch')
    V.assign_properties('EALL', CPS3, mat, t=.001)

    # ASSIGN HOMOGENEOUS BCS TO MODEL
    V.assign_prescribed_bc(10, (X,Y))

    # CREATE LOAD STAGE AND TO IT ASSIGN INHOMOGENEOUS BCS
    stage = V.create_static_stage(solver=NEWTON)
    stage.assign_prescribed_bc(20, X, .24e-3)
    stage.assign_prescribed_bc(20, Y, .12e-3)
    stage.assign_prescribed_bc(30, X,  .3e-3)
    stage.assign_prescribed_bc(30, Y, .24e-3)
    stage.assign_prescribed_bc(40, X, .06e-3)
    stage.assign_prescribed_bc(40, Y, .12e-3)

    # RUN THE STAGE
    stage.run()

    # CHECK THE RESULTS AGAINST ANALYTIC SOLUTION
    stage = V.stages.last
    field = stage.increments[-1].field_outputs['S']
    error = []
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.333333333), error.append('Sxx')
        assert allclose(data[:,1], 1333.333333333), error.append('Syy')
        assert allclose(data[:,2],  400.), error.append('Sxy')

    field = stage.increments[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3), error.append('Exx')
        assert allclose(data[:,1], 1e-3), error.append('Eyy')
        assert allclose(data[:,2], 1e-3), error.append('Exy')

    if error:
        error = ','.join(error)
        message = 'PATCH TEST FAILED, FAILED RESULTS: {0}'.format(error)
        logging.error(message)
        raise Exception(message)
    else:
        logging.info('PATCH TEST PASSED')

    if plot:
        # VISUALIZE RESULTS
        V.Plot2D(show=1, deformed=1)

    # CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
    V = fe_model(mesh=mesh, jobid='Tri3PlaneStressPatch')
    V.assign_properties('EALL', CPS3, mat, t=.001)

    # ASSIGN HOMOGENEOUS BCS TO MODEL
    V.assign_prescribed_bc(10, (X,Y))

    # CREATE LOAD STAGE AND TO IT ASSIGN INHOMOGENEOUS BCS
    stage = V.create_static_stage()
    stage.assign_prescribed_bc(20, X, .24e-3)
    stage.assign_prescribed_bc(20, Y, .12e-3)
    stage.assign_prescribed_bc(30, X,  .3e-3)
    stage.assign_prescribed_bc(30, Y, .24e-3)
    stage.assign_prescribed_bc(40, X, .06e-3)
    stage.assign_prescribed_bc(40, Y, .12e-3)

    # RUN THE STAGE
    stage.run()

    # CHECK THE RESULTS AGAINST ANALYTIC SOLUTION
    stage = V.stages.last
    field = stage.increments[-1].field_outputs['S']
    error = []
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1333.333333333), error.append('Sxx')
        assert allclose(data[:,1], 1333.333333333), error.append('Syy')
        assert allclose(data[:,2],  400.), error.append('Sxy')

    field = stage.increments[-1].field_outputs['E']
    for value in field.values:
        data = value.data
        assert allclose(data[:,0], 1e-3), error.append('Exx')
        assert allclose(data[:,1], 1e-3), error.append('Eyy')
        assert allclose(data[:,2], 1e-3), error.append('Exy')

    if error:
        error = ','.join(error)
        message = 'PATCH TEST FAILED, FAILED RESULTS: {0}'.format(error)
        logging.error(message)
        raise Exception(message)
    else:
        logging.info('PATCH TEST PASSED')

    if plot:
        # VISUALIZE RESULTS
        V.Plot2D(show=1, deformed=1)

def test(plot=False):
    PlaneStressTria3PatchDemo(plot=plot)

if __name__ == '__main__':
    test(plot=True)
