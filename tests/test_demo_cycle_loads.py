from felab import *

def CycleLoads(plot=False):
    mesh = unit_square_mesh(nx=1, ny=1)
    mat = Material('Mat-1', elastic={'E': 500, 'Nu': 0})

    V = fe_model(mesh=mesh)
    V.create_element_block('Block-1', ALL)
    V.assign_properties('Block-1', CPS4, mat)

    V.fix_nodes(ILO)

    step = V.create_static_step(solver=NEWTON)
    step.assign_prescribed_bc(IHI, X, .1)
    step.run()
    if plot:
        V.Plot2D(show=1, deformed=True)

    step = V.create_static_step(solver=NEWTON)
    step.assign_prescribed_bc(IHI, X, 0)
    step.run()
    if plot:
        V.Plot2D(show=1, deformed=True)

    step = V.create_static_step(solver=NEWTON)
    step.remove_bc(IHI, X)
    step.assign_surface_load(IHI, [100, 0])
    step.run()
    if plot:
        V.Plot2D(show=1, deformed=True)


def test(plot=False):
    CycleLoads(plot=plot)


if __name__ == '__main__':
    test(plot=True)
