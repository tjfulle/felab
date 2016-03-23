from pyfem2 import *

# READ MESH
mesh = AbaqusMesh(filename='EC3SFP1.inp')

# CREATE MATERIAL MODEL
mat = Material('Material-1', elastic={'E':1e6, 'Nu':.25})

# CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
V = FiniteElementModel(mesh=mesh, jobid='PlaneStressTria3Patch')
V.AssignProperties('EALL', PlaneStressTria3, mat, t=.001)

# ASSIGN HOMOGENEOUS BCS TO MODEL
V.PrescribedBC(10, (X,Y))

# CREATE LOAD STEP AND TO IT ASSIGN INHOMOGENEOUS BCS
step = V.StaticStep(solver=NEWTON)
step.PrescribedBC(20, X, .24e-3)
step.PrescribedBC(20, Y, .12e-3)
step.PrescribedBC(30, X,  .3e-3)
step.PrescribedBC(30, Y, .24e-3)
step.PrescribedBC(40, X, .06e-3)
step.PrescribedBC(40, Y, .12e-3)

# RUN THE STEP
step.run()

# CHECK THE RESULTS AGAINST ANALYTIC SOLUTION
step = V.steps.last
field = step.frames[-1].field_outputs['S']
error = []
for value in field.values:
    data = value.data
    assert allclose(data[:,0], 1333.333333333), error.append('Sxx')
    assert allclose(data[:,1], 1333.333333333), error.append('Syy')
    assert allclose(data[:,2],  400.), error.append('Sxy')

field = step.frames[-1].field_outputs['E']
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

if not os.getenv('NOGRAPHICS'):
    # VISUALIZE RESULTS
    V.Plot2D(show=1, deformed=1)

# CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
V = FiniteElementModel(mesh=mesh, jobid='PlaneStressTria3Patch')
V.AssignProperties('EALL', PlaneStressTria3, mat, t=.001)

# ASSIGN HOMOGENEOUS BCS TO MODEL
V.PrescribedBC(10, (X,Y))

# CREATE LOAD STEP AND TO IT ASSIGN INHOMOGENEOUS BCS
step = V.StaticStep()
step.PrescribedBC(20, X, .24e-3)
step.PrescribedBC(20, Y, .12e-3)
step.PrescribedBC(30, X,  .3e-3)
step.PrescribedBC(30, Y, .24e-3)
step.PrescribedBC(40, X, .06e-3)
step.PrescribedBC(40, Y, .12e-3)

# RUN THE STEP
step.run()

# CHECK THE RESULTS AGAINST ANALYTIC SOLUTION
step = V.steps.last
field = step.frames[-1].field_outputs['S']
error = []
for value in field.values:
    data = value.data
    assert allclose(data[:,0], 1333.333333333), error.append('Sxx')
    assert allclose(data[:,1], 1333.333333333), error.append('Syy')
    assert allclose(data[:,2],  400.), error.append('Sxy')

field = step.frames[-1].field_outputs['E']
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

if not os.getenv('NOGRAPHICS'):
    # VISUALIZE RESULTS
    V.Plot2D(show=1, deformed=1)
