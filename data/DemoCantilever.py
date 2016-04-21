from pyfem2 import *

def Cantilever8():
    filename = 'Cantilever8.inp'
    mesh = AbaqusMesh(filename)
    mat = Material('Material-1', elastic={'E':100., 'Nu':.3})
    V = FiniteElementModel(mesh=mesh)
    V.AssignProperties('EAll', PlaneStressQuad8, mat, t=1)
    V.FixNodes((1, 22, 33))
    step = V.StaticStep()
    step.ConcentratedLoad(49, Y, .01)
    step.run(increments=10)
    V.WriteResults()
Cantilever8()
