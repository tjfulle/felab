import pytest
from numpy import *
from conf import *
from felab import *


def test_print_matrix_1():
    eps = lambda: random.random()
    eps = lambda: 0
    p = array(
        [
            [0 + eps(), 2 + eps()],
            [1 + eps(), 2 + eps()],
            [2 + eps(), 2 + eps()],
            [0 + eps(), 1 + eps()],
            [1 + eps(), 1 + eps()],
            [2 + eps(), 1 + eps()],
            [0 + eps(), 0 + eps()],
            [1 + eps(), 0 + eps()],
            [2 + eps(), 0 + eps()],
        ],
        dtype=float64,
    )
    t = array(
        [
            [0, 4, 1],
            [1, 5, 2],
            [1, 4, 5],
            [0, 3, 4],
            [3, 7, 4],
            [4, 8, 5],
            [4, 7, 8],
            [3, 6, 7],
        ],
        dtype=int,
    )

    k, h, Too = 12, 250, 25
    V = fe_model(jobid="test-1")
    V.mesh = Mesh(p=p, t=t)
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(k)
    V.create_element_block("ElementBlock1", ALL)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    V.assign_initial_temperature(ALL, 50)
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(BOUNDARY, T, 50)
    fun = lambda x: ones_like(x[:, 0])
    stage.HeatGeneration(ALL, fun)
    stage.run()
    stage.print_stiffness_structure(style="numeric", stream=sys.stdout)
    V.write_results()


if __name__ == "__main__":
    test_print_matrix_1()
