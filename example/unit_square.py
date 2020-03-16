from felab.fe_model import fe_model
from felab.elemlib import CPS4
from felab.constants import ALL, ILO, IHI, X
from felab.material import Material
from felab.mesh import unit_square_mesh


def demo_unit_square(plot=False):
    mesh = unit_square_mesh(nx=2, ny=2)
    mat = Material("Mat", elastic={"E": 1000, "Nu": 0})

    V = fe_model(mesh=mesh)
    V.element_block("All", ALL)
    V.assign_properties("All", CPS4, mat)
    V.fix_nodes(ILO)

    step = V.static_step()
    step.assign_prescribed_bc(IHI, X, 0.1)
    step.run()
    step.print_stiffness_structure(style="latex")
    if plot:
        V.Plot2D(deformed=1, show=1)

    step = V.static_step()
    step.assign_prescribed_bc(IHI, X, 0)
    step.run()

    if plot:
        V.Plot2D(deformed=1, show=1)


if __name__ == "__main__":
    demo_unit_square(plot=True)
