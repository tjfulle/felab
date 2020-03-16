from felab.fe_model import fe_model
from felab.elemlib import CPE4
from felab.constants import ILO, ALL


def test_gravity_load(plot=False):
    V = fe_model(jobid="Gravity")
    V.rectilinear_mesh(nx=20, ny=5, lx=100, ly=10)
    V.create_material("Material-1")
    V.materials["Material-1"].Density(1.0)
    V.materials["Material-1"].elastic(E=10e6, Nu=0.333)
    V.create_element_block("Block1", ALL)
    V.assign_properties("Block1", CPE4, "Material-1", t=1)
    stage = V.create_static_stage()
    stage.fix_nodes(ILO)
    stage.assign_gravity_load(ALL, [0, -1e2])
    stage.run()
    V.write_results()

    if plot:
        ax = V.Plot2D(label="Undeformed", color="orange")
        V.Plot2D(deformed=1, label="Deformed", color="blue", ax=ax, show=1)


if __name__ == "__main__":
    test_gravity_load(plot=True)
