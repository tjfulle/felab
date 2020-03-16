from numpy import nan
from felab.fe_model import fe_model
from felab.elemlib import B2D2, L2D2
from felab.material import Material
from felab.constants import X, Y, TZ
from felab.mesh import Mesh


def test_beam_column():
    nodtab = [[1, -4, 3], [2, 0, 0], [3, 0, 3], [4, nan, nan], [5, 4, 3]]
    eletab = [[1, 1, 3], [2, 3, 5], [3, 1, 2], [4, 2, 3], [5, 2, 5]]
    mesh = Mesh(nodtab=nodtab, eletab=eletab)

    Ec, Em = 30000, 200000
    mat1 = Material("Material-1", elastic={"E": Ec, "Nu": 0.3})
    mat2 = Material("Material-2", elastic={"E": Em, "Nu": 0.3})

    V = fe_model(mesh=mesh, jobid="B2D2")
    V.create_element_block("B1", (1, 2))
    V.create_element_block("B2", (3, 5))
    V.create_element_block("B3", (4,))
    V.assign_properties("B1", B2D2, mat1, A=0.02, Izz=0.004)
    V.assign_properties("B2", L2D2, mat2, A=0.001)
    V.assign_properties("B3", L2D2, mat2, A=0.003)

    V.assign_prescribed_bc(1, (X, Y, TZ))
    V.assign_prescribed_bc(5, Y)

    step = V.create_static_step()
    step.assign_concentrated_load(2, Y, 100)
    step.assign_concentrated_load(5, TZ, 200)
    step.assign_concentrated_load(5, X, 300)
    step.run()
    V.write_results()


if __name__ == "__main__":
    demo_beam_column()
