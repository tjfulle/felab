from pyfem2 import *

def Cantilever8():
    filename = 'Cantilever8.inp'
    mesh = AbaqusMesh(filename)
    V = FiniteElementModel(mesh=mesh)
    for eb in mesh.element_blocks.values():
        V.AssignProperties(eb.name, eb.eletyp, eb.material, t=1)
    V.FixNodes((1, 22, 33))
    step = V.StaticStep(solver=NEWTON, increments=10)
    step.ConcentratedLoad(49, Y, .01)
    step.run()
    V.WriteResults()
Cantilever8()
