import os
import pytest
import numpy as np

try:
    import distmesh as dm
except ImportError:
    dm = None
from felab.fe_model import fe_model
from felab.elemlib import DC2D3
from felab.constants import T, ALL, ILO, IHI, JLO, JHI, BOUNDARY


def test_heat_transfer_1(data_path):
    def solution(x, N=20):
        def fun(k):
            a = np.sin(k * np.pi * (1 + x[:, 0]) / 2.0) / (k ** 3 * np.sinh(k * np.pi))
            b = np.sinh(k * np.pi * (1 + x[:, 1]) / 2.0) + np.sinh(
                k * np.pi * (1 - x[:, 1]) / 2.0
            )
            return a * b

        u = (1 - x[:, 0] ** 2) / 2.0 - 16.0 / np.pi ** 3 * sum(
            [fun(k) for k in range(1, N, 2)], 0
        )
        return u

    V = fe_model()
    V.genesis_mesh(os.path.join(data_path, "UniformPlateTria3Fine.g"))
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(1.0)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    step = V.create_heat_transfer_step()
    step.HeatGeneration(ALL, 1)
    step.assign_prescribed_bc(BOUNDARY, T, 0)
    step.run()
    Tn = solution(V.mesh.coord)
    err = np.sqrt(np.mean((step.dofs.flatten() - Tn) ** 2)) / np.sqrt(np.mean(Tn ** 2))
    assert err < 1e-4


def test_heat_transfer_2(data_path):
    def solution(x):
        return 2.0 * (1.0 + x[:, 1]) / ((3.0 + x[:, 0]) ** 2 + (1 + x[:, 1]) ** 2)

    V = fe_model()
    V.genesis_mesh(os.path.join(data_path, "UniformPlateTria3.g"))
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(1.0)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(BOUNDARY, T, solution)
    step.HeatGeneration(ALL, 0)
    step.run()
    Tn = solution(V.mesh.coord)
    err = np.sqrt(np.mean((step.dofs.flatten() - Tn) ** 2)) / np.sqrt(np.mean(Tn ** 2))
    assert err < 5e-5


@pytest.mark.skipif(dm is None, reason="distmesh not import")
def test_heat_transfer_3():
    if dm is None:
        return
    def solution(x, N=20):
        def fun(n):
            return (
                np.sin(2.0 * n * np.pi / 3.0)
                / n ** 2
                / np.sinh(n * np.pi)
                * np.sin(n * np.pi * x[:, 0])
                * np.sinh(n * np.pi * x[:, 1])
            )

        return 450.0 / np.pi ** 2 * sum([fun(n) for n in range(1, N)], 0)

    np.random.seed(190)  # Always the same results
    fd = lambda p: dm.drectangle(p, 0, 1, 0, 1)
    fh = dm.huniform
    coord, elecon = dm.distmesh2d(
        fd, fh, 0.05, (0, 0, 1, 1), [(0, 0), (0, 1), (1, 0), (1, 1)]
    )
    f2 = lambda x: np.where(x[:, 0] <= 2.0 / 3.0, 75 * x[:, 0], 150 * (1 - x[:, 0]))
    V = fe_model()
    V.create_mesh(p=coord, t=elecon)
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(1.0)
    V.create_element_block("ElementBlock1", ALL)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(JLO, T, 0)
    step.assign_prescribed_bc(JHI, T, f2)
    step.assign_prescribed_bc(ILO, T, 0)
    step.assign_prescribed_bc(IHI, T, 0)
    step.HeatGeneration(ALL, 0)
    step.run()
    Tn = solution(V.mesh.coord)
    err = np.sqrt(np.mean((step.dofs.flatten() - Tn) ** 2)) / np.sqrt(np.mean(Tn ** 2))
    assert err < 5e-3


def test_heat_transfer_4():
    def solution(x, q0=1.0, k=1.0, N=20):
        def fun(n):
            al = 0.5 * (2.0 * n - 1.0) * np.pi
            top = (-1) ** n * np.cos(al * x[:, 1]) * np.cosh(al * x[:, 0])
            bot = al ** 3 * np.cosh(al)
            return top / bot

        return (
            q0
            / 2
            / k
            * ((1 - x[:, 1] ** 2) + 4.0 * sum([fun(n) for n in range(1, N)], 0))
        )

    nodtab = [
        [1, 0.0, 0.0],
        [2, 0.5, 0.0],
        [3, 1.0, 0.0],
        [4, 0.0, 0.5],
        [5, 0.5, 0.5],
        [6, 1.0, 0.5],
        [7, 0.0, 1.0],
        [8, 0.5, 1.0],
        [9, 1.0, 1.0],
    ]
    eletab = [
        [1, 1, 2, 5],
        [2, 1, 5, 4],
        [3, 2, 3, 6],
        [4, 2, 6, 5],
        [5, 4, 5, 8],
        [6, 4, 8, 7],
        [7, 5, 6, 9],
        [8, 5, 9, 8],
    ]
    V = fe_model()
    V.create_mesh(nodtab=nodtab, eletab=eletab)
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(1.0)
    V.create_element_block("ElementBlock1", ALL)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(IHI, T)
    step.assign_prescribed_bc(JHI, T)
    step.HeatGeneration(ALL, 1)
    step.run()
    Tn = solution(V.mesh.coord)
    err = np.sqrt(np.mean((step.dofs.flatten() - Tn) ** 2)) / np.sqrt(np.mean(Tn ** 2))
    assert err < 0.045


@pytest.mark.skipif(dm is None, reason="distmesh not import")
def test_heat_transfer_5():
    if dm is None:
        return
    def solution(x, q0=1.0, k=1.0, N=20):
        def fun(n):
            al = 0.5 * (2.0 * n - 1.0) * np.pi
            top = (-1) ** n * np.cos(al * x[:, 1]) * np.cosh(al * x[:, 0])
            bot = al ** 3 * np.cosh(al)
            return top / bot

        return (
            q0
            / 2
            / k
            * ((1 - x[:, 1] ** 2) + 4.0 * sum([fun(n) for n in range(1, N)], 0))
        )

    np.random.seed(190)  # Always the same results
    fd = lambda p: dm.drectangle(p, 0, 1, 0, 1)
    fh = dm.huniform
    coord, elecon = dm.distmesh2d(
        fd, fh, 0.025, (0, 0, 1, 1), [(0, 0), (0, 1), (1, 0), (1, 1)]
    )
    V = fe_model()
    V.create_mesh(p=coord, t=elecon)
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(1.0)
    V.create_element_block("ElementBlock1", ALL)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(IHI, T)
    step.assign_prescribed_bc(JHI, T)
    step.HeatGeneration(ALL, 1)
    step.run()
    Tn = solution(V.mesh.coord)
    err = np.sqrt(np.mean((step.dofs.flatten() - Tn) ** 2)) / np.sqrt(np.mean(Tn ** 2))
    assert err < 1e-4
