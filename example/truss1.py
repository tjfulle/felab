from felab import *


def demo_truss():
    # Create the mesh from tables of nodes and elements
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

    # Create the model problem
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    V = FEModel(jobid="Truss1", mesh=mesh)

    # Create a material and define the elastic properties
    mat = Material(name="Material-1", elastic=dict(E=1000, Nu=0.29))

    # Define an element block of 3D 2-node link elements
    V.element_block(name="ElementBlock1", elements=ALL)

    # Assign material properties
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
    el = Element(type="L2D2")
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, A=A
    )

    # Apply boundary conditions
    V.dirichlet_bc(1, (X, Y))
    V.dirichlet_bc(12, Y)

    # Apply concentrated loads
    step = V.static_step()
    step.concentrated_load((3, 5, 9, 11), Y, -10)
    step.concentrated_load(7, Y, -16)

    step.run()
    V.write_results()


if __name__ == "__main__":
    demo_truss()
