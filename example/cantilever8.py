from felab.fe_model import FEModel
from felab.mesh import abaqus_mesh
from felab.material import Material
from felab.elemlib import CPS8
from felab.constants import Y


def demo_cantilever8():
    filename = "./data/Cantilever8.inp"
    mesh = abaqus_mesh(filename)
    mat = Material("Material-1", elastic={"E": 100.0, "Nu": 0.3})
    V = FEModel(mesh=mesh)
    V.assign_properties(element_block="EAll", element_type=CPS8, material=mat, t=1)
    V.fix_nodes((1, 22, 33))
    step = V.static_step()
    step.concentrated_load(49, Y, 0.01)
    step.run(frames=10)
    V.write_results()


if __name__ == "__main__":
    demo_cantilever8()
