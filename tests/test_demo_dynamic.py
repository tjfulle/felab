from felab import *

def DynamicLoadStage(plot=False):
    mesh = unit_square_mesh(nx=1, ny=1)
    mat = Material('Mat-1', elastic={'E': 500, 'Nu': 0})
    mat.Density(1.)

    V = fe_model(mesh=mesh)
    V.create_element_block('Block-1', ALL)
    V.assign_properties('Block-1', CPE4, mat)

    stage = V.create_dynamic_stage(period=1e-6, increments=10)
    stage.assign_prescribed_bc(IHI, X, .1)
    stage.run()
    if plot:
        V.Plot2D(show=1, deformed=True)


def test(plot=False):
    DynamicLoadStage(plot=plot)


if __name__ == '__main__':
    test(plot=True)
