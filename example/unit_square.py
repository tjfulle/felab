from felab import *
from felab.elemlib import CPS4
from felab.io.plot import plot2d


def demo_unit_square(plot=False):
    mesh = unit_square_mesh(nx=2, ny=2)
    mat = Material(name="Mat", elastic=dict(E=1000, Nu=0))

    V = FEModel(mesh=mesh)
    V.element_block(name="All", elements=ALL)
    V.assign_properties(element_block="All", element_type=CPS4, material=mat)
    V.fix_nodes(ILO)

    step = V.static_step()
    step.dirichlet_bc(IHI, X, 0.1)
    step.run()
    step.print_stiffness_structure(style="latex")
    if plot:
        plot2d(model=V, deformed=1, show=1)

    step = V.static_step()
    step.dirichlet_bc(IHI, X, 0)
    step.run()

    if plot:
        plot2d(model=V, deformed=1, show=1)


if __name__ == "__main__":
    demo_unit_square(plot=True)
