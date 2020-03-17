from felab.fe_model import FEModel
from felab.elemlib import CPE4
from felab.constants import ILO, ALL
from felab.io.plot import plot2d


def test_gravity_load(plot=False):
    V = FEModel(jobid="Gravity")
    V.rectilinear_mesh(nx=20, ny=5, lx=100, ly=10)
    V.material("Material-1")
    V.materials["Material-1"].Density(1.0)
    V.materials["Material-1"].elastic(E=10e6, Nu=0.333)
    V.element_block(name="Block1", elements=ALL)
    V.assign_properties(
        element_block="Block1", element_type=CPE4, material="Material-1", t=1
    )
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
