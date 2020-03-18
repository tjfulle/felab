from felab import *


def demo_cantilever8():
    filename = "./data/Cantilever8.inp"
    mesh = abaqus_mesh(filename)
    mat = Material(name="Material-1", elastic=dict(E=100.0, Nu=0.3))
    V = FEModel(mesh=mesh)
    el = Element(type="CPS8")
    V.assign_properties(element_block="EAll", element_type=el, material=mat, t=1)
    V.fix_nodes((1, 22, 33))
    step = V.static_step()
    step.concentrated_load(49, Y, 0.01)
    step.run(frames=10)
    V.write_results()


if __name__ == "__main__":
    demo_cantilever8()
