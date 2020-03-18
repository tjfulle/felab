from felab import *
from felab.elemlib import CPE4
from felab.io.plot import plot2d


def demo_beam_quad(plot=False):
    mesh = rectilinear_mesh2d(nx=10, ny=2, lx=10, ly=2)
    mat = Material(name="Material-1", elastic=dict(E=20000, Nu=0))

    V = FEModel(mesh=mesh, jobid="QuadBeam")
    V.element_block(name="ElementBlock1", elements=ALL)
    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE4, material=mat, t=1
    )

    step = V.static_step()
    step.fix_dofs(ILO)
    step = V.static_step()
    step.concentrated_load(IHI, Y, -10)
    step.run()
    V.write_results()
    if plot:
        plot2d(model=V, show=1, deformed=1)


if __name__ == "__main__":
    demo_beam_quad(plot=True)
