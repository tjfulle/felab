from felab.fe_model import fe_model
from felab.elemlib import CPE4
from felab.constants import ALL, X, IHI
from felab.mesh import unit_square_mesh
from felab.material import Material
from felab.io.plot import plot2d


def demo_dynamic_load_step(plot=False):
    mesh = unit_square_mesh(nx=1, ny=1)
    mat = Material("Mat-1", elastic={"E": 500, "Nu": 0})
    mat.Density(1.0)

    V = fe_model(mesh=mesh)
    V.element_block(name="Block-1", elements=ALL)
    V.assign_properties(element_block="Block-1", element_type=CPE4, material=mat)

    step = V.dynamic_step(period=1e-6, frames=10)
    step.assign_prescribed_bc(IHI, X, 0.1)
    step.run()
    if plot:
        plot2d(model=V, show=1, deformed=True)


if __name__ == "__main__":
    demo_dynamic_load_step(plot=True)
