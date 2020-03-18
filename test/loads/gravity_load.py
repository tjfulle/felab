from felab import *
from felab.io.plot import plot2d


def test_gravity_load(plot=False):
    mesh = rectilinear_mesh(nx=20, ny=5, lx=100, ly=10)
    V = FEModel(jobid="Gravity", mesh=mesh)
    el = Element(type="CPE4")
    mat = Material(name="Material-1", density=1.0, elastic=dict(E=10e6, Nu=0.333))
    V.element_block(name="Block1", elements=ALL)
    V.assign_properties(element_block="Block1", element_type=el, material=mat, t=1)
    step = V.static_step()
    step.fix_nodes(ILO)
    step.gravity_load(ALL, [0, -1e2])
    step.run()
    V.write_results()

    if plot:
        ax = plot2d(model=V, label="Undeformed", color="orange")
        plot2d(model=V, deformed=1, label="Deformed", color="blue", ax=ax, show=1)


if __name__ == "__main__":
    test_gravity_load(plot=True)
