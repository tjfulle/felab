from pyfem2 import *
V = Plane2DModel(jobid='PlaneStrainQuad4Patch')
V.AbaqusMesh(filename='EC4SFP1.inp')
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
#V.Plot2D(show=1)
V.Solve()
V.WriteResults()

# Average stress must be 1600 in x and y
step = V.steps.last
field = step.frames[-1].field_outputs['S']
for value in field.values:
    data = value.data
    assert allclose(data[:,0], 1600.), 'Wrong Sxx'
    assert allclose(data[:,1], 1600.), 'Wrong Syy'
    assert allclose(data[:,2],  800.), 'Wrong Szz'
    assert allclose(data[:,3],  400.), 'Wrong Sxy'
