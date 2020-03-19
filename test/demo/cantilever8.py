from felab import *


def demo_cantilever8(data_path):
    filename = os.path.join(data_path, "Cantilever8.inp")
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
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_cantilever8(data_path)
