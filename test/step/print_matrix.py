import sys
import numpy as np
from felab import *
from felab.elemlib import DC2D3


def test_print_matrix_1():
    eps = lambda: np.random.random()
    eps = lambda: 0
    p = np.array(
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
        dtype=float,
    )
    t = np.array(
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
    mesh = Mesh(p=p, t=t)
    V = FEModel(jobid="test-1", mesh=mesh)

    mat = Material(name="Material-1", isotropic_thermal_conductivity=k)
    V.element_block(name="ElementBlock1", elements=ALL)
    V.assign_properties(element_block="ElementBlock1", element_type=DC2D3, material=mat)
    V.initial_temperature(ALL, 50)
    step = V.heat_transfer_step()
    step.dirichlet_bc(BOUNDARY, T, 50)
    fun = lambda x: np.ones_like(x[:, 0])
    step.HeatGeneration(ALL, fun)
    step.run()
    step.print_stiffness_structure(style="numeric", stream=sys.stdout)
    V.write_results()


if __name__ == "__main__":
    test_print_matrix_1()
