import numpy as np
import felab.util.tty as tty
from felab.fe_model import fe_model
from felab.elemlib import CPS3
from felab.constants import NEWTON, X, Y
from felab.mesh import abaqus_mesh
from felab.material import Material


def demo_plane_stress_tria3_patch(plot=False):
    # READ MESH
    mesh = abaqus_mesh(filename="./data/EC3SFP1.inp")

    # CREATE MATERIAL MODEL
    mat = Material("Material-1", elastic={"E": 1e6, "Nu": 0.25})

    # CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
    V = fe_model(mesh=mesh, jobid="Tri3PlaneStressPatch")
    V.assign_properties("EALL", CPS3, mat, t=0.001)

    # ASSIGN HOMOGENEOUS BCS TO MODEL
    V.assign_prescribed_bc(10, (X, Y))

    # CREATE LOAD STEP AND TO IT ASSIGN INHOMOGENEOUS BCS
    step = V.static_step(solver=NEWTON)
    step.assign_prescribed_bc(20, X, 0.24e-3)
    step.assign_prescribed_bc(20, Y, 0.12e-3)
    step.assign_prescribed_bc(30, X, 0.3e-3)
    step.assign_prescribed_bc(30, Y, 0.24e-3)
    step.assign_prescribed_bc(40, X, 0.06e-3)
    step.assign_prescribed_bc(40, Y, 0.12e-3)

    # RUN THE STEP
    step.run()

    # CHECK THE RESULTS AGAINST ANALYTIC SOLUTION
    step = V.steps.last
    field = step.frames[-1].field_outputs["S"]
    error = []
    for value in field.values:
        data = value.data
        assert np.allclose(data[:, 0], 1333.333333333), error.append("Sxx")
        assert np.allclose(data[:, 1], 1333.333333333), error.append("Syy")
        assert np.allclose(data[:, 2], 400.0), error.append("Sxy")

    field = step.frames[-1].field_outputs["E"]
    for value in field.values:
        data = value.data
        assert np.allclose(data[:, 0], 1e-3), error.append("Exx")
        assert np.allclose(data[:, 1], 1e-3), error.append("Eyy")
        assert np.allclose(data[:, 2], 1e-3), error.append("Exy")

    if error:
        error = ",".join(error)
        message = "PATCH TEST FAILED, FAILED RESULTS: {0}".format(error)
        tty.die(message)
    else:
        tty.info("PATCH TEST PASSED")

    if plot:
        # VISUALIZE RESULTS
        V.Plot2D(show=1, deformed=1)

    # CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
    V = fe_model(mesh=mesh, jobid="Tri3PlaneStressPatch")
    V.assign_properties("EALL", CPS3, mat, t=0.001)

    # ASSIGN HOMOGENEOUS BCS TO MODEL
    V.assign_prescribed_bc(10, (X, Y))

    # CREATE LOAD STEP AND TO IT ASSIGN INHOMOGENEOUS BCS
    step = V.static_step()
    step.assign_prescribed_bc(20, X, 0.24e-3)
    step.assign_prescribed_bc(20, Y, 0.12e-3)
    step.assign_prescribed_bc(30, X, 0.3e-3)
    step.assign_prescribed_bc(30, Y, 0.24e-3)
    step.assign_prescribed_bc(40, X, 0.06e-3)
    step.assign_prescribed_bc(40, Y, 0.12e-3)

    # RUN THE STEP
    step.run()

    # CHECK THE RESULTS AGAINST ANALYTIC SOLUTION
    step = V.steps.last
    field = step.frames[-1].field_outputs["S"]
    error = []
    for value in field.values:
        data = value.data
        assert np.allclose(data[:, 0], 1333.333333333), error.append("Sxx")
        assert np.allclose(data[:, 1], 1333.333333333), error.append("Syy")
        assert np.allclose(data[:, 2], 400.0), error.append("Sxy")

    field = step.frames[-1].field_outputs["E"]
    for value in field.values:
        data = value.data
        assert np.allclose(data[:, 0], 1e-3), error.append("Exx")
        assert np.allclose(data[:, 1], 1e-3), error.append("Eyy")
        assert np.allclose(data[:, 2], 1e-3), error.append("Exy")

    if error:
        error = ",".join(error)
        message = "PATCH TEST FAILED, FAILED RESULTS: {0}".format(error)
        tty.die(message)
    else:
        tty.info("PATCH TEST PASSED")

    if plot:
        # VISUALIZE RESULTS
        V.Plot2D(show=1, deformed=1)


if __name__ == "__main__":
    demo_plane_stress_tria3_patch(plot=True)
