from felab import *
from felab.io.plot import plot2d


def demo_plane_elasticity2(data_path, plot=False):

    # READ MESH FROM FILE
    mesh = genesis_mesh(os.path.join(data_path, "PlateWithHoleTria3.g"))

    # ELASTIC MATERIAL MODEL
    mat = Material(name="Material-1", elastic=dict(E=10e6, Nu=0.29))

    # CREATE THE MODEL AND ASSIGN PROPERTIES
    el = Element(type="CPS3")
    V = FEModel(mesh=mesh, jobid="Plane2")
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
    )

    # PRESCRIBE FIXED BCS TO MODEL
    V.dirichlet_bc("LeftHandSide", X)
    V.dirichlet_bc("BottomLeft", Y)

    # CREATE LOAD STEP AND PRESCRIBED NONHOMOGENEOUS BCS TO IT
    step = V.static_step()
    step.dirichlet_bc("RightHandSide", X, 0.1)

    # RUN THE STEP TO SOLVE FOR THE UNKNOWN DEGREES OF FREEDOM
    step.run()

    # WRITE RESULTS
    V.write_results()

    if plot:
        # VISUALIZE RESULTS
        plot2d(model=V, show=1, deformed=1)


if __name__ == "__main__":
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_plane_elasticity2(data_path, plot=True)
