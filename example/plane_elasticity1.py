#!/usr/bin/env python
from felab.fe_model import fe_model
from felab.elemlib import CPE4
from felab.constants import X, Y
from felab.io.plot import plot2d


def demo_plane_elasticity1(plot=False):
    V = fe_model(jobid="Plane1")
    V.genesis_mesh("./data/PlateWithHoleQuad4.g")

    mat = V.material("Material-1")
    mat.elastic(E=10e6, Nu=0.29)

    V.assign_properties("ElementBlock1", CPE4, mat, t=1)

    step = V.static_step()
    step.assign_prescribed_bc("LeftHandSide", X, 0.0)
    step.assign_prescribed_bc("PinNode", Y, 0.0)
    step.assign_prescribed_bc("RightHandSide", X, 0.1)

    step.run()

    V.write_results()

    if plot:
        plot2d(model=V, show=1, deformed=1, colorby="Ux")


if __name__ == "__main__":
    demo_plane_elasticity1(plot=True)
