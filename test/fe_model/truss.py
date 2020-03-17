from numpy import allclose, array, nan
from felab.fe_model import FEModel
from felab.elemlib import L3D2, L2D2, B2D2
from felab.constants import X, Y, Z, ALL, TZ


def test_model_truss_0():
    nodtab = [[1, 0, 0], [2, 3, 4], [3, 0, 4]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=70e9, Nu=0.3)
    V.element_block(name="B1", elements=ALL)
    V.assign_properties(
        element_block="B1", element_type=L2D2, material="Material-1", A=5 * 0.01 * 0.01
    )
    step = V.static_step()
    step.dirichlet_bc(1, X, -0.05)
    step.dirichlet_bc((2, 3), (X, Y))
    step.concentrated_load(1, Y, 1000e3)
    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data
    assert allclose([[-0.05, 0.0882842], [0.0, 0.0], [0.0, 0.0]], u)


def test_model_truss_1():
    nodtab = [
        [1, 0, 0, 0],
        [2, 10, 5, 0],
        [3, 10, 0, 0],
        [4, 20, 8, 0],
        [5, 20, 0, 0],
        [6, 30, 9, 0],
        [7, 30, 0, 0],
        [8, 40, 8, 0],
        [9, 40, 0, 0],
        [10, 50, 5, 0],
        [11, 50, 0, 0],
        [12, 60, 0, 0],
    ]
    eletab = [
        [1, 1, 3],
        [2, 3, 5],
        [3, 5, 7],
        [4, 7, 9],
        [5, 9, 11],
        [6, 11, 12],
        [7, 1, 2],
        [8, 2, 4],
        [9, 4, 6],
        [10, 6, 8],
        [11, 8, 10],
        [12, 10, 12],
        [13, 2, 3],
        [14, 4, 5],
        [15, 6, 7],
        [16, 8, 9],
        [17, 10, 11],
        [18, 2, 5],
        [19, 4, 7],
        [20, 7, 8],
        [21, 9, 10],
    ]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=1000, Nu=0.333)
    Abot, Atop, Abat, Adia = 2, 10, 3, 1
    A = [
        Abot,
        Abot,
        Abot,
        Abot,
        Abot,
        Abot,
        Atop,
        Atop,
        Atop,
        Atop,
        Atop,
        Atop,
        Abat,
        Abat,
        Abat,
        Abat,
        Abat,
        Adia,
        Adia,
        Adia,
        Adia,
    ]
    V.element_block(name="element_block-1", elements=ALL)
    V.assign_properties(
        element_block="element_block-1", element_type=L3D2, material="Material-1", A=A
    )

    step = V.static_step()
    step.dirichlet_bc(1, (X, Y))
    step.dirichlet_bc(12, Y)
    step.dirichlet_bc(ALL, Z)

    step.concentrated_load((3, 5, 9, 11), Y, -10)
    step.concentrated_load(7, Y, -16)

    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data
    RF = V.steps.last.frames[-1].field_outputs["RF"].data
    assert allclose(
        [
            [0.0, 0.0, 0.0],
            [0.80954, -1.7756, 0.0],
            [0.28, -1.79226, 0.0],
            [0.899, -2.29193, 0.0],
            [0.56, -2.3166, 0.0],
            [0.8475, -2.38594, 0.0],
            [0.8475, -2.42194, 0.0],
            [0.796, -2.29193, 0.0],
            [1.135, -2.3166, 0.0],
            [0.88546, -1.7756, 0.0],
            [1.415, -1.79226, 0.0],
            [1.695, 0.0, 0.0],
        ],
        u,
    )
    assert allclose(
        [
            [0.0, 28.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [-0.0, -0.0, 0.0],
            [-0.0, 0.0, 0.0],
            [-0.0, 0.0, 0.0],
            [0.0, -0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 28.0, 0.0],
        ],
        RF,
    )


def test_model_truss_2():
    nodtab = [
        [1, 0, 0],
        [2, 10, 5],
        [3, 10, 0],
        [4, 20, 8],
        [5, 20, 0],
        [6, 30, 9],
        [7, 30, 0],
        [8, 40, 8],
        [9, 40, 0],
        [10, 50, 5],
        [11, 50, 0],
        [12, 60, 0],
    ]
    eletab = [
        [1, 1, 3],
        [2, 3, 5],
        [3, 5, 7],
        [4, 7, 9],
        [5, 9, 11],
        [6, 11, 12],
        [7, 1, 2],
        [8, 2, 4],
        [9, 4, 6],
        [10, 6, 8],
        [11, 8, 10],
        [12, 10, 12],
        [13, 2, 3],
        [14, 4, 5],
        [15, 6, 7],
        [16, 8, 9],
        [17, 10, 11],
        [18, 2, 5],
        [19, 4, 7],
        [20, 7, 8],
        [21, 9, 10],
    ]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=1000, Nu=0.333)
    Abot, Atop, Abat, Adia = 2, 10, 3, 1
    A = [
        Abot,
        Abot,
        Abot,
        Abot,
        Abot,
        Abot,
        Atop,
        Atop,
        Atop,
        Atop,
        Atop,
        Atop,
        Abat,
        Abat,
        Abat,
        Abat,
        Abat,
        Adia,
        Adia,
        Adia,
        Adia,
    ]
    V.element_block(name="element_block-1", elements=ALL)
    V.assign_properties(
        element_block="element_block-1", element_type=L2D2, material="Material-1", A=A
    )

    step = V.static_step()
    step.dirichlet_bc(1, (X, Y))
    step.dirichlet_bc(12, Y)

    step.concentrated_load((3, 5, 9, 11), Y, -10)
    step.concentrated_load(7, Y, -16)

    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data
    RF = V.steps.last.frames[-1].field_outputs["RF"].data
    assert allclose(
        [
            [0.0, 0.0],
            [0.80954, -1.7756],
            [0.28, -1.79226],
            [0.899, -2.29193],
            [0.56, -2.3166],
            [0.8475, -2.38594],
            [0.8475, -2.42194],
            [0.796, -2.29193],
            [1.135, -2.3166],
            [0.88546, -1.7756],
            [1.415, -1.79226],
            [1.695, 0.0],
        ],
        u,
    )
    assert allclose(
        [
            [0.0, 28.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [-0.0, -0.0],
            [-0.0, 0.0],
            [-0.0, 0.0],
            [0.0, -0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 28.0],
        ],
        RF,
    )


def test_model_truss_3():
    nodtab = [[1, 0.0, 0.0], [2, 3.0, 4.0], [3, 0.0, 4.0]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=70e9, Nu=0.333)
    A = 5 * 0.01 * 0.01
    V.element_block(name="element_block-1", elements=ALL)
    V.assign_properties(
        element_block="element_block-1", element_type=L2D2, material="Material-1", A=A
    )
    step = V.static_step()
    step.fix_nodes((2, 3))
    step.dirichlet_bc(1, X, -0.05)
    step.concentrated_load(1, Y, 1000e3)
    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data
    assert allclose([[-0.05, 0.0882842], [0.0, 0.0], [0.0, 0.0]], u)


def test_model_truss_4():
    nodtab = [[1, 0.0, 0.0, 0.0], [2, 3.0, 4.0, 0.0], [3, 0.0, 4.0, 0.0]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=70e9, Nu=0.333)
    A = 5 * 0.01 * 0.01
    V.element_block(name="element_block-1", elements=ALL)
    V.assign_properties(
        element_block="element_block-1", element_type=L3D2, material="Material-1", A=A
    )
    step = V.static_step()
    step.fix_nodes((2, 3))
    step.dirichlet_bc(1, X, -0.05)
    step.dirichlet_bc(1, Z)
    step.concentrated_load(1, Y, 1000e3)
    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data
    assert allclose([[-0.05, 0.0882842, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], u)


def test_model_truss_5():
    nodtab = [[1, 72, 0, 0], [2, 0, 36, 0], [3, 0, 36, 72], [4, 0, 0, -48]]
    eletab = [[1, 1, 2], [2, 1, 3], [3, 1, 4]]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=10e4, Nu=0.333)
    A = [0.302, 0.729, 0.187]
    V.element_block(name="element_block-1", elements=ALL)
    V.assign_properties(
        element_block="element_block-1", element_type=L3D2, material="Material-1", A=A
    )
    step = V.static_step()
    # Boundary conditions
    step.fix_nodes((2, 3, 4))
    step.dirichlet_bc(1, Y)
    # Concentrated force in 'z' direction on node 1
    step.concentrated_load(1, Z, -1000)
    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data
    ua = (
        array(
            [
                [-8.5337228, 0.0, -31.9486913],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        )
        / 10
    )
    assert allclose(u, ua)


def test_model_truss_6():
    # Set up problem space
    nodtab = [
        [1, -37.5, 0, 200],
        [2, 37.5, 0, 200],
        [3, -37.5, 37.5, 100],
        [4, 37.5, 37.5, 100],
        [5, 37.5, -37.5, 100],
        [6, -37.5, -37.5, 100],
        [7, -100, 100, 0],
        [8, 100, 100, 0],
        [9, 100, -100, 0],
        [10, -100, -100, 0],
    ]
    eletab = [
        [1, 1, 2],
        [2, 1, 4],
        [3, 2, 3],
        [4, 1, 5],
        [5, 2, 6],
        [6, 2, 4],
        [7, 2, 5],
        [8, 1, 3],
        [9, 1, 6],
        [10, 3, 6],
        [11, 4, 5],
        [12, 3, 4],
        [13, 5, 6],
        [14, 3, 10],
        [15, 6, 7],
        [16, 4, 9],
        [17, 5, 8],
        [18, 4, 7],
        [19, 3, 8],
        [20, 5, 10],
        [21, 6, 9],
        [22, 6, 10],
        [23, 3, 7],
        [24, 5, 9],
        [25, 4, 8],
    ]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)

    # Define element blocks
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=10e6, Nu=0.333)
    A = [
        0.033,
        2.015,
        2.015,
        2.015,
        2.015,
        2.823,
        2.823,
        2.823,
        2.823,
        0.01,
        0.01,
        0.014,
        0.014,
        0.98,
        0.98,
        0.98,
        0.98,
        1.76,
        1.76,
        1.76,
        1.76,
        2.44,
        2.44,
        2.44,
        2.44,
    ]
    V.element_block(name="element_block-1", elements=ALL)
    V.assign_properties(
        element_block="element_block-1", element_type=L3D2, material="Material-1", A=A
    )

    step = V.static_step()

    # Define boundary conditons
    step.fix_nodes([7, 8, 9, 10])

    # Define concentrated loads
    step.concentrated_load(1, X, 1000)
    step.concentrated_load(1, Y, 10000)
    step.concentrated_load(1, Z, -5000)
    step.concentrated_load(2, Y, 10000)
    step.concentrated_load(2, Z, -5000)
    step.concentrated_load(3, X, 500)
    step.concentrated_load(6, X, 500)

    # Solve and write results
    step.run()
    u = V.steps.last.frames[-1].field_outputs["U"].data

    assert allclose(
        [
            [0.00851510679597, 0.349956039184, -0.0221277138856],
            [0.0319156311642, 0.349956039184, -0.0322420125936],
            [0.0115296378344, -0.00976991195147, -0.108526000994],
            [-0.00403948591251, -0.00878106640481, -0.115392665916],
            [0.000447704186587, -0.00508916981698, 0.0705078974296],
            [0.00704244773528, -0.00410032427032, 0.0773745623519],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
        u,
    )


def test_model_truss_beam_0():
    nodtab = [[1, -4, 3], [2, 0, 0], [3, 0, 3], [4, nan, nan], [5, 4, 3]]
    eletab = [[1, 1, 3], [2, 3, 5], [3, 1, 2], [4, 2, 3], [5, 2, 5]]
    V = FEModel()
    V.ne_mesh(nodtab=nodtab, eletab=eletab)
    Ec, Em = 30000, 200000
    V.material("Material-1")
    V.materials["Material-1"].elastic(E=Ec, Nu=0.3)
    V.material("Material-2")
    V.materials["Material-2"].elastic(E=Em, Nu=0.3)
    V.element_block(name="B1", elements=(1, 2))
    V.element_block(name="B2", elements=(3, 5))
    V.element_block(name="B3", elements=(4,))
    V.assign_properties(
        element_block="B1", element_type=B2D2, material="Material-1", A=0.02, Izz=0.004
    )
    V.assign_properties(
        element_block="B2", element_type=L2D2, material="Material-2", A=0.001
    )
    V.assign_properties(
        element_block="B3", element_type=L2D2, material="Material-2", A=0.003
    )
    V.dirichlet_bc(1, (X, Y, TZ))
    V.dirichlet_bc(5, Y)
    step = V.static_step()
    step.concentrated_load(2, Y, 100)
    step.concentrated_load(5, TZ, 200)
    step.concentrated_load(5, X, 300)
    step.run()
