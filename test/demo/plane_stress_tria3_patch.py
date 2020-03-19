import numpy as np
from felab import *
import felab.util.tty as tty
from felab.io.plot import plot2d


def demo_plane_stress_tria3_patch(data_path, plot=False):
    # READ MESH
    mesh = abaqus_mesh(filename=os.path.join(data_path, "EC3SFP1.inp"))

    # CREATE MATERIAL MODEL
    mat = Material(name="Material-1", elastic=dict(E=1e6, Nu=0.25))

    # CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
    el = Element(type="CPS3")
    V = FEModel(mesh=mesh, jobid="Tri3PlaneStressPatch")
    V.assign_properties(element_block="EALL", element_type=el, material=mat, t=0.001)

    # ASSIGN HOMOGENEOUS BCS TO MODEL
    V.dirichlet_bc(10, (X, Y))

    # CREATE LOAD STEP AND TO IT ASSIGN INHOMOGENEOUS BCS
    step = V.static_step(solver=NEWTON)
    step.dirichlet_bc(20, X, 0.24e-3)
    step.dirichlet_bc(20, Y, 0.12e-3)
    step.dirichlet_bc(30, X, 0.3e-3)
    step.dirichlet_bc(30, Y, 0.24e-3)
    step.dirichlet_bc(40, X, 0.06e-3)
    step.dirichlet_bc(40, Y, 0.12e-3)

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
        plot2d(model=V, show=1, deformed=1)

    # CREATE FINITE ELEMENT MODEL AND ASSIGN PROPERTIES
    el = Element(type="CPS3")
    V = FEModel(mesh=mesh, jobid="Tri3PlaneStressPatch")
    V.assign_properties(element_block="EALL", element_type=el, material=mat, t=0.001)

    # ASSIGN HOMOGENEOUS BCS TO MODEL
    V.dirichlet_bc(10, (X, Y))

    # CREATE LOAD STEP AND TO IT ASSIGN INHOMOGENEOUS BCS
    step = V.static_step()
    step.dirichlet_bc(20, X, 0.24e-3)
    step.dirichlet_bc(20, Y, 0.12e-3)
    step.dirichlet_bc(30, X, 0.3e-3)
    step.dirichlet_bc(30, Y, 0.24e-3)
    step.dirichlet_bc(40, X, 0.06e-3)
    step.dirichlet_bc(40, Y, 0.12e-3)

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
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_plane_stress_tria3_patch(data_path, plot=True)
