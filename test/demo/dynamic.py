from felab import *
from felab.io.plot import plot2d


def demo_dynamic_load_step(plot=False):
    mesh = unit_square_mesh(nx=1, ny=1)
    mat = Material(name="Mat-1", density=1.0, elastic=dict(E=500, Nu=0))

    V = FEModel(mesh=mesh)
    el = Element(type="CPE4")
    V.element_block(name="Block-1", elements=ALL)
    V.assign_properties(element_block="Block-1", element_type=el, material=mat)

    step = V.dynamic_step(period=1e-6, frames=10)
    step.dirichlet_bc(IHI, X, 0.1)
    step.run()
    if plot:
        plot2d(model=V, show=1, deformed=True)


if __name__ == "__main__":
    demo_dynamic_load_step(plot=True)
