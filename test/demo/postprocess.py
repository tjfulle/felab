from felab import *
from felab.io.exodusii import EXOFileReader


def demo_postprocess(data_path):
    mesh = genesis_mesh(os.path.join(data_path, "PlateWithHoleQuad4.g"))
    mat = Material(name="Material-1", elastic=dict(E=100, Nu=0.2))

    V = FEModel(mesh=mesh, jobid="PlateWithHoleQuad4")

    el = Element(type="CPE4")
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
    )
    V.dirichlet_bc("LeftHandSide", X)
    V.fix_nodes("PinNode")

    step = V.static_step()
    step.surface_load(IHI, [1, 0])
    step.run()

    V.write_results()

    F = EXOFileReader("PlateWithHoleQuad4.exo")
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


if __name__ == "__main__":
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_postprocess(data_path)
