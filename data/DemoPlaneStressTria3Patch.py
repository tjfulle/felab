from pyfem2 import *
V = FiniteElementModel(jobid='PlaneStressTria3Patch')
V.AbaqusMesh(filename='EC3SFP1.inp')
mat = V.Material('Material-1')
mat.Elastic(E=1e6, Nu=.25)
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
logging.info('PATCH TEST PASSED')
