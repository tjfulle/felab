import numpy as np
import felab.util.tty as tty
from felab.fe_model import fe_model
from felab.elemlib import CPE4
from felab.constants import X, Y


def demo_plane_stran_quad4_patch(plot=False):
    V = fe_model(jobid="Quad4PlaneStrainPatch")
    V.abaqus_mesh(filename="./data/EC4SFP1.inp")
    mat = V.create_material("Material-1")
    mat.elastic(E=1e6, Nu=0.25)
    V.assign_properties("EALL", CPE4, mat, t=0.001)

    step = V.create_static_step()
    step.assign_prescribed_bc(10, (X, Y), 0.0)
    step.assign_prescribed_bc(20, X, 0.24e-3)
    step.assign_prescribed_bc(20, Y, 0.12e-3)
    step.assign_prescribed_bc(30, X, 0.3e-3)
    step.assign_prescribed_bc(30, Y, 0.24e-3)
    step.assign_prescribed_bc(40, X, 0.06e-3)
    step.assign_prescribed_bc(40, Y, 0.12e-3)
    step.run()
    if plot:
        V.Plot2D(show=1)
    V.write_results()

    # Average stress must be 1600 in x and y
    step = V.steps.last
    field = step.frames[-1].field_outputs["S"]
    for value in field.values:
        data = value.data
        assert np.allclose(data[:, 0], 1600.0), "Wrong Sxx"
        assert np.allclose(data[:, 1], 1600.0), "Wrong Syy"
        assert np.allclose(data[:, 2], 800.0), "Wrong Szz"
        assert np.allclose(data[:, 3], 400.0), "Wrong Sxy"
    tty.info("PATCH TEST PASSED")


if __name__ == "__main__":
    demo_plane_stran_quad4_patch(plot=True)
