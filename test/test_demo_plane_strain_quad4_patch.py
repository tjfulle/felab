import logging
from felab import *


def PlaneStrainQuad4PatchDemo(plot=False):
    V = fe_model(jobid="Quad4PlaneStrainPatch")
    V.abaqus_mesh(filename="./data/EC4SFP1.inp")
    mat = V.create_material("Material-1")
    mat.elastic(E=1e6, Nu=0.25)
    V.assign_properties("EALL", CPE4, mat, t=0.001)

    stage = V.create_static_stage()
    stage.assign_prescribed_bc(10, (X, Y), 0.0)
    stage.assign_prescribed_bc(20, X, 0.24e-3)
    stage.assign_prescribed_bc(20, Y, 0.12e-3)
    stage.assign_prescribed_bc(30, X, 0.3e-3)
    stage.assign_prescribed_bc(30, Y, 0.24e-3)
    stage.assign_prescribed_bc(40, X, 0.06e-3)
    stage.assign_prescribed_bc(40, Y, 0.12e-3)
    stage.run()
    if plot:
        V.Plot2D(show=1)
    V.write_results()

    # Average stress must be 1600 in x and y
    stage = V.stages.last
    field = stage.increments[-1].field_outputs["S"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1600.0), "Wrong Sxx"
        assert allclose(data[:, 1], 1600.0), "Wrong Syy"
        assert allclose(data[:, 2], 800.0), "Wrong Szz"
        assert allclose(data[:, 3], 400.0), "Wrong Sxy"
    logging.info("PATCH TEST PASSED")


def test(plot=False):
    PlaneStrainQuad4PatchDemo(plot=plot)


if __name__ == "__main__":
    test(plot=True)
