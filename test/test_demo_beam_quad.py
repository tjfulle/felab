from felab import *


def BeamQuad(plot=False):
    mesh = rectilinear_mesh2d(nx=10, ny=2, lx=10, ly=2)
    mat = Material("Material-1", elastic={"E": 20000, "Nu": 0})

    V = fe_model(mesh=mesh, jobid="QuadBeam")
    V.create_element_block("ElementBlock1", ALL)
    V.assign_properties("ElementBlock1", CPE4, mat, t=1)

    stage = V.create_static_stage()
    stage.fix_dofs(ILO)
    stage = V.create_static_stage()
    stage.assign_concentrated_load(IHI, Y, -10)
    stage.run()
    V.write_results()
    if plot:
        V.Plot2D(show=1, deformed=1)


def test(plot=False):
    BeamQuad(plot=plot)


if __name__ == "__main__":
    test(plot=True)
