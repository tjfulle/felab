#!/usr/bin/env python
from numpy import *
from felab import *

def PlaneElasticity1(plot=False):
    V = fe_model(jobid='Plane1')
    V.genesis_mesh('./data/PlateWithHoleQuad4.g')

    mat = V.create_material('Material-1')
    mat.elastic(E=10e6, Nu=.29)

    V.assign_properties('ElementBlock1', CPE4, mat, t=1)

    step = V.create_static_step()
    step.assign_prescribed_bc('LeftHandSide', X, 0.)
    step.assign_prescribed_bc('PinNode', Y, 0.)
    step.assign_prescribed_bc('RightHandSide', X, .1)

    step.run()

    V.write_results()

    if plot:
        V.Plot2D(show=1, deformed=1, colorby='Ux')

def test(plot=False):
    PlaneElasticity1(plot=plot)

if __name__ == '__main__':
    test(plot=True)
