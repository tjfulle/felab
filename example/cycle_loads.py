from felab import *
from felab.io.plot import plot2d


def demo_cycle_loads(plot=False):
    mesh = unit_square_mesh(nx=1, ny=1)
    mat = Material(name="Mat-1", elastic=dict(E=500, Nu=0))
    el = Element(type="CPS4")

    V = FEModel(mesh=mesh)
    V.element_block(name="Block-1", elements=ALL)
    V.assign_properties(element_block="Block-1", element_type=el, material=mat)

    V.fix_nodes(ILO)

    step = V.static_step(solver=NEWTON)
    step.dirichlet_bc(IHI, X, 0.1)
    step.run()
    if plot:
        plot2d(model=V, show=1, deformed=True)

    step = V.static_step(solver=NEWTON)
    step.dirichlet_bc(IHI, X, 0)
    step.run()
    if plot:
        plot2d(model=V, show=1, deformed=True)

    step = V.static_step(solver=NEWTON)
    step.remove_bc(IHI, X)
    step.surface_load(IHI, [100, 0])
    step.run()
    if plot:
        plot2d(model=V, show=1, deformed=True)


if __name__ == "__main__":
    demo_cycle_loads(plot=True)
