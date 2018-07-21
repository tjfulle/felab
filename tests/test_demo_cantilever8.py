from felab import *

def Cantilever8():
    filename = './data/Cantilever8.inp'
    mesh = abaqus_mesh(filename)
    mat = Material('Material-1', elastic={'E':100., 'Nu':.3})
    V = fe_model(mesh=mesh)
    V.assign_properties('EAll', CPS8, mat, t=1)
    V.fix_nodes((1, 22, 33))
    step = V.create_static_step()
    step.assign_concentrated_load(49, Y, .01)
    step.run(increments=10)
    V.write_results()


def test():
    Cantilever8()


if __name__ == '__main__':
    test()
