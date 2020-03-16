#!/usr/bin/env python
from felab.fe_model import fe_model
from felab.elemlib import CPE4
from felab.constants import X, Y


def test_plane_elasticity1(plot=False):
    V = fe_model(jobid="Plane1")
    V.genesis_mesh("./data/PlateWithHoleQuad4.g")

    mat = V.create_material("Material-1")
    mat.elastic(E=10e6, Nu=0.29)

    V.assign_properties("ElementBlock1", CPE4, mat, t=1)

    stage = V.create_static_stage()
    stage.assign_prescribed_bc("LeftHandSide", X, 0.0)
    stage.assign_prescribed_bc("PinNode", Y, 0.0)
    stage.assign_prescribed_bc("RightHandSide", X, 0.1)

    stage.run()

    V.write_results()

    if plot:
        V.Plot2D(show=1, deformed=1, colorby="Ux")


if __name__ == "__main__":
    demo_plane_elasticity1(plot=True)
