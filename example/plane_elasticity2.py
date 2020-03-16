from felab.fe_model import fe_model
from felab.mesh import genesis_mesh
from felab.elemlib import CPS3
from felab.material import Material
from felab.constants import X, Y


def test_plane_elasticity2(plot=False):

    # READ MESH FROM FILE
    mesh = genesis_mesh("./data/PlateWithHoleTria3.g")

    # ELASTIC MATERIAL MODEL
    mat = Material("Material-1", elastic={"E": 10e6, "Nu": 0.29})

    # CREATE THE MODEL AND ASSIGN PROPERTIES
    V = fe_model(mesh=mesh, jobid="Plane2")
    V.assign_properties("ElementBlock1", CPS3, mat, t=1)

    # PRESCRIBE FIXED BCS TO MODEL
    V.assign_prescribed_bc("LeftHandSide", X)
    V.assign_prescribed_bc("BottomLeft", Y)

    # CREATE LOAD STEP AND PRESCRIBED NONHOMOGENEOUS BCS TO IT
    step = V.create_static_step()
    step.assign_prescribed_bc("RightHandSide", X, 0.1)

    # RUN THE STEP TO SOLVE FOR THE UNKNOWN DEGREES OF FREEDOM
    step.run()

    # WRITE RESULTS
    V.write_results()

    if plot:
        # VISUALIZE RESULTS
        V.Plot2D(show=1, deformed=1)


if __name__ == "__main__":
    demo_plane_elasticity2(plot=True)
