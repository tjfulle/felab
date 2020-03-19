#!/usr/bin/env python
from felab import *


def demo_truss():

    # Create mesh and define the finite element model
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
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    V = FEModel(mesh=mesh, jobid="Truss2")

    # Define element blocks
    mat = Material(name="Material-1", elastic=dict(E=10e6, Nu=0.333))
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
    el = Element(type="L3D2")
    V.element_block(name="ElementBlock1", elements=ALL)
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, A=A
    )

    # Define boundary conditons
    V.fix_nodes((7, 8, 9, 10))

    # Define concentrated loads
    step = V.static_step()

    P1, P2, P3, P4 = 1000, 10000, -5000, 500
    step.concentrated_load(1, X, P1)
    step.concentrated_load(1, Y, P2)
    step.concentrated_load(1, Z, P3)
    step.concentrated_load(2, Y, P2)
    step.concentrated_load(2, Z, P3)
    step.concentrated_load((3, 6), X, P4)

    # Solve and write results
    step.run()
    V.write_results()


if __name__ == "__main__":
    demo_truss()
