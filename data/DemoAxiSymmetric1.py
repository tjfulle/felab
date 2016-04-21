#!/usr/bin/env python
from pyfem2 import *
import matplotlib.pyplot as plt

def AxisymmetricRing():
    V = FiniteElementModel()
    nodtab = [[1, 1000.0, 0.0],
              [2, 1002.0, 0.0],
              [3, 1002.0, 1.0],
              [4, 1000.0, 1.0]]
    eletab = [[1, 1, 2, 3, 4]]
    V.Mesh(nodtab=nodtab, eletab=eletab)

    # MATERIAL
    mat = V.Material('Material-1')
    mat.Elastic(E=30e6, Nu=.3)

    V.ElementBlock('EALL', (1,))
    V.AssignProperties('EALL', AxiSymmetricQuad4, mat, formulation=1)
    V.PrescribedBC(1, Zr)
    V.PrescribedBC(2, Zr)

    step = V.StaticStep()
    step.Pressure((1, S2), 1000)
    step.Pressure((1, S3), 1000)
    step.Pressure((1, S4), 1000)
    step.run()

    s = step.frames[-1].field_outputs['S']
    #print(s.data[:,:,[0,1,2]])
    #print(s.data[:,:,3])
    assert allclose(s.data[:,:,[0,1,2]], -1000.)
    assert allclose(s.data[:,:,3], 0)

    e = step.frames[-1].field_outputs['E']
    assert allclose(e.data[:,:,[0,1,2]], -1.3333e-5)
    assert allclose(e.data[:,:,3], 0)

AxisymmetricRing()
