from felab.fe_model import fe_model
from felab.constants import NEWTON, ALL, X, ILO, IHI
from felab.material import Material
from felab.elemlib import CPS4
from felab.mesh import unit_square_mesh


def test_cycle_loads(plot=False):
    mesh = unit_square_mesh(nx=1, ny=1)
    mat = Material("Mat-1", elastic={"E": 500, "Nu": 0})

    V = fe_model(mesh=mesh)
    V.create_element_block("Block-1", ALL)
    V.assign_properties("Block-1", CPS4, mat)

    V.fix_nodes(ILO)

    stage = V.create_static_stage(solver=NEWTON)
    stage.assign_prescribed_bc(IHI, X, 0.1)
    stage.run()
    if plot:
        V.Plot2D(show=1, deformed=True)

    stage = V.create_static_stage(solver=NEWTON)
    stage.assign_prescribed_bc(IHI, X, 0)
    stage.run()
    if plot:
        V.Plot2D(show=1, deformed=True)

    stage = V.create_static_stage(solver=NEWTON)
    stage.remove_bc(IHI, X)
    stage.assign_surface_load(IHI, [100, 0])
    stage.run()
    if plot:
        V.Plot2D(show=1, deformed=True)


if __name__ == "__main__":
    demo_cycle_loads(plot=True)
