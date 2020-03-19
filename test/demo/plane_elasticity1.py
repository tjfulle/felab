from felab import *
from felab.io.plot import plot2d


def demo_plane_elasticity1(data_path, plot=False):
    mesh = genesis_mesh(os.path.join(data_path, "PlateWithHoleQuad4.g"))
    V = FEModel(jobid="Plane1", mesh=mesh)

    el = Element(type="CPE4")
    mat = Material(name="Material-1", elastic=dict(E=10e6, Nu=0.29))
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
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
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_plane_elasticity1(data_path, plot=True)
