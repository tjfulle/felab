from felab import *
from felab.elemlib import B2D2, L2D2


def demo_beam_column():
    nodtab = [[1, -4, 3], [2, 0, 0], [3, 0, 3], [4, None, None], [5, 4, 3]]
    eletab = [[1, 1, 3], [2, 3, 5], [3, 1, 2], [4, 2, 3], [5, 2, 5]]
    mesh = Mesh(nodtab=nodtab, eletab=eletab)

    Ec, Em = 30000, 200000
    mat1 = Material(name="Material-1", elastic=dict(E=Ec, Nu=0.3))
    mat2 = Material(name="Material-2", elastic=dict(E=Em, Nu=0.3))

    V = FEModel(mesh=mesh, jobid="B2D2")
    V.element_block(name="B1", elements=(1, 2))
    V.element_block(name="B2", elements=(3, 5))
    V.element_block(name="B3", elements=(4,))
    el1 = Element(type="B2D2")
    V.assign_properties(
        element_block="B1", element_type=el1, material=mat1, A=0.02, Izz=0.004
    )
    el2 = Element(type="L2D2")
    V.assign_properties(element_block="B2", element_type=el2, material=mat2, A=0.001)
    V.assign_properties(element_block="B3", element_type=el2, material=mat2, A=0.003)

    V.dirichlet_bc(1, (X, Y, TZ))
    V.dirichlet_bc(5, Y)

    step = V.static_step()
    step.concentrated_load(2, Y, 100)
    step.concentrated_load(5, TZ, 200)
    step.concentrated_load(5, X, 300)
    step.run()
    V.write_results()


if __name__ == "__main__":
    demo_beam_column()
