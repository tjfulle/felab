from felab import *
from felab.elemlib import CPE4
from felab.io.exodusii import EXOFileReader
from felab.io.plot import plot2d
from felab.material import Material


def demo_quarter_plate(plot=False):
    mesh = genesis_mesh("./data/PlateWithHoleQuad4QuarterSym.g")
    V = FEModel(jobid="PlateWithHoleQuad4QuarterSym", mesh=mesh)
    mat = Material(name="Material-1", elastic=dict(E=100, Nu=0.2))
    V.assign_properties(element_block="", element_type=CPE4, material=mat, t=1)
    V.dirichlet_bc("SymYZ", X)
    V.dirichlet_bc("SymXZ", Y)
    V.initial_temperature(ALL, 60)

    step = V.static_step("Step-1")
    step.surface_load("RightHandSide", [1, 0])
    step.run()

    V.write_results()

    F = EXOFileReader("PlateWithHoleQuad4QuarterSym.exo")
    max_p = [0.0, None]
    max_u = [0.0, None]
    for step in F.steps.values():
        for frame in step.frames:
            u = frame.field_outputs["U"]
            for value in u.values:
                u1 = value.magnitude
                if max_u[0] < u1:
                    max_u = [u1, value]

            s = frame.field_outputs["S"]
            for value in s.values:
                s1 = value.max_principal
                if max_p[0] < s1:
                    max_p = [s1, value]

    # External and internal element numbers
    if plot:
        plot2d(model=V, deformed=1, show=True)


if __name__ == "__main__":
    demo_quarter_plate(plot=True)
