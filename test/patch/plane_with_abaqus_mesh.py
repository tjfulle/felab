import os
from numpy import allclose
from felab.constants import X, Y
from felab.fe_model import fe_model
from felab.material import Material
from felab.elemlib import CPE4, CPS4, CPS3, CPE3


def test_quad4_plane_strain(data_path):
    V = fe_model()
    V.abaqus_mesh(filename=os.path.join(data_path, "EC4SFP1.inp"))
    mat = Material("Material-1", elastic={"E": 1e6, "Nu": 0.25})
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

    # Average stress must be 1600 in x and y
    field = step.increments[-1].field_outputs["S"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1600.0)
        assert allclose(data[:, 1], 1600.0)
        assert allclose(data[:, 2], 800.0)
        assert allclose(data[:, 3], 400.0)
    field = step.increments[-1].field_outputs["E"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1e-3)
        assert allclose(data[:, 1], 1e-3)
        assert allclose(data[:, 3], 1e-3)


def test_quad4_plane_stress(data_path):
    V = fe_model()
    V.abaqus_mesh(filename=os.path.join(data_path, "EC4SFP1.inp"))
    mat = Material("Material-1", elastic={"E": 1e6, "Nu": 0.25})
    V.assign_properties("EALL", CPS4, mat, t=0.001)
    step = V.create_static_step()
    step.assign_prescribed_bc(10, (X, Y), 0.0)
    step.assign_prescribed_bc(20, X, 0.24e-3)
    step.assign_prescribed_bc(20, Y, 0.12e-3)
    step.assign_prescribed_bc(30, X, 0.3e-3)
    step.assign_prescribed_bc(30, Y, 0.24e-3)
    step.assign_prescribed_bc(40, X, 0.06e-3)
    step.assign_prescribed_bc(40, Y, 0.12e-3)
    step.run()
    field = step.increments[-1].field_outputs["S"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1333.33333333)
        assert allclose(data[:, 1], 1333.33333333)
        assert allclose(data[:, 2], 400.0)
    field = step.increments[-1].field_outputs["E"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1e-3)
        assert allclose(data[:, 1], 1e-3)
        assert allclose(data[:, 2], 1e-3)


def test_tria3_plane_stress(data_path):
    V = fe_model(jobid="PlaneStressTria3Patch")
    V.abaqus_mesh(filename=os.path.join(data_path, "EC3SFP1.inp"))
    mat = Material("Material-1", elastic={"E": 1e6, "Nu": 0.25})
    V.assign_properties("EALL", CPS3, mat, t=0.001)
    step = V.create_static_step()
    step.assign_prescribed_bc(10, (X, Y), 0.0)
    step.assign_prescribed_bc(20, X, 0.24e-3)
    step.assign_prescribed_bc(20, Y, 0.12e-3)
    step.assign_prescribed_bc(30, X, 0.3e-3)
    step.assign_prescribed_bc(30, Y, 0.24e-3)
    step.assign_prescribed_bc(40, X, 0.06e-3)
    step.assign_prescribed_bc(40, Y, 0.12e-3)
    step.run()
    V.write_results()
    field = step.increments[-1].field_outputs["S"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1333.333333333), "Wrong Sxx"
        assert allclose(data[:, 1], 1333.333333333), "Wrong Syy"
        assert allclose(data[:, 2], 400.0), "Wrong Sxy"
    field = step.increments[-1].field_outputs["E"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1e-3)
        assert allclose(data[:, 1], 1e-3)
        assert allclose(data[:, 2], 1e-3)


def test_tria3_plane_strain(data_path):
    V = fe_model()
    V.abaqus_mesh(filename=os.path.join(data_path, "EC3SFP1.inp"))
    mat = Material("Material-1", elastic={"E": 1e6, "Nu": 0.25})
    V.assign_properties("EALL", CPE3, mat, t=0.001)
    step = V.create_static_step()
    step.assign_prescribed_bc(10, (X, Y), 0.0)
    step.assign_prescribed_bc(20, X, 0.24e-3)
    step.assign_prescribed_bc(20, Y, 0.12e-3)
    step.assign_prescribed_bc(30, X, 0.3e-3)
    step.assign_prescribed_bc(30, Y, 0.24e-3)
    step.assign_prescribed_bc(40, X, 0.06e-3)
    step.assign_prescribed_bc(40, Y, 0.12e-3)
    step.run()
    # Average stress must be 1600 in x and y
    field = step.increments[-1].field_outputs["S"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1600.0)
        assert allclose(data[:, 1], 1600.0)
        assert allclose(data[:, 2], 800.0)
        assert allclose(data[:, 3], 400.0)
    field = step.increments[-1].field_outputs["E"]
    for value in field.values:
        data = value.data
        assert allclose(data[:, 0], 1e-3)
        assert allclose(data[:, 1], 1e-3)
        assert allclose(data[:, 3], 1e-3)
