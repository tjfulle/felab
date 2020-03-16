from felab.fe_model import fe_model
from felab.mesh import abaqus_mesh
from felab.material import Material
from felab.elemlib import CPS8
from felab.constants import Y


def test_cantilever8():
    filename = "./data/Cantilever8.inp"
    mesh = abaqus_mesh(filename)
    mat = Material("Material-1", elastic={"E": 100.0, "Nu": 0.3})
    V = fe_model(mesh=mesh)
    V.assign_properties("EAll", CPS8, mat, t=1)
    V.fix_nodes((1, 22, 33))
    step = V.create_static_step()
    step.assign_concentrated_load(49, Y, 0.01)
    step.run(increments=10)
    V.write_results()


if __name__ == "__main__":
    demo_cantilever8()
