from numpy import allclose
from felab.fe_model import FEModel
from felab.elemlib import CAX4
from felab.constants import Zr, S2, S3, S4


def demo_axisymmetric_ring():
    V = FEModel()
    nodtab = [[1, 1000.0, 0.0], [2, 1002.0, 0.0], [3, 1002.0, 1.0], [4, 1000.0, 1.0]]
    eletab = [[1, 1, 2, 3, 4]]
    V.ne_mesh(nodtab=nodtab, eletab=eletab)

    # MATERIAL
    mat = V.material("Material-1")
    mat.elastic(E=30e6, Nu=0.3)

    V.element_block(name="EALL", elements=(1,))
    V.assign_properties(
        element_block="EALL", element_type=CAX4, material=mat, formulation=1
    )
    V.dirichlet_bc(1, Zr)
    V.dirichlet_bc(2, Zr)

    step = V.static_step()
    step.pressure((1, S2), 1000)
    step.pressure((1, S3), 1000)
    step.pressure((1, S4), 1000)
    step.run()

    s = step.frames[-1].field_outputs["S"]
    # print(s.data[:,:,[0,1,2]])
    # print(s.data[:,:,3])
    assert allclose(s.data[:, :, [0, 1, 2]], -1000.0)
    assert allclose(s.data[:, :, 3], 0)

    e = step.frames[-1].field_outputs["E"]
    assert allclose(e.data[:, :, [0, 1, 2]], -1.3333e-5)
    assert allclose(e.data[:, :, 3], 0)


if __name__ == "__main__":
    demo_axisymmetric_ring()
