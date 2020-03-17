#!/usr/bin/env python
from felab.fe_model import FEModel
from felab.elemlib import CPE4
from felab.constants import X, Y
from felab.io.plot import plot2d


def demo_plane_elasticity1(plot=False):
    V = FEModel(jobid="Plane1")
    V.genesis_mesh("./data/PlateWithHoleQuad4.g")

    mat = V.material("Material-1")
    mat.elastic(E=10e6, Nu=0.29)

    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE4, material=mat, t=1
    )

    step = V.static_step()
    step.dirichlet_bc("LeftHandSide", X, 0.0)
    step.dirichlet_bc("PinNode", Y, 0.0)
    step.dirichlet_bc("RightHandSide", X, 0.1)

    step.run()

    V.write_results()

    if plot:
        plot2d(model=V, show=1, deformed=1, colorby="Ux")


if __name__ == "__main__":
    demo_plane_elasticity1(plot=True)
