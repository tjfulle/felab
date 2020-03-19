import numpy as np
from felab import *
from felab.io.plot import plot2d


def demo_volume_locking(data_path, plot=False):
    mu = 1.0
    for (i, Nu) in enumerate((0.0, 0.2, 0.45, 0.499)):
        E = 2.0 * mu * (1.0 + Nu)

        # Analytic solution
        mesh = Mesh(filename=os.path.join(data_path, "QuarterCylinderQuad4.g"))
        p = 1.0
        a, b = mesh.coord[0, 1], mesh.coord[-1, 0]
        u = np.zeros_like(mesh.coord)
        for (j, x) in enumerate(mesh.coord):
            r = np.sqrt(x[0] ** 2 + x[1] ** 2)
            term1 = (1.0 + Nu) * a ** 2 * b ** 2 * p / (E * (b ** 2 - a ** 2))
            term2 = 1.0 / r + (1.0 - 2.0 * Nu) * r / b ** 2
            ur = term1 * term2
            u[j, :] = ur * x[:] / r
        mesh.put_nodal_solution("VolumeLocking.Analytic", u)
        if plot:
            ax = plot2d(
                mesh=mesh,
                xy=u + mesh.coord,
                color="orange",
                label=r"Analytic, $\nu=%g$" % Nu,
                weight=8,
            )

        # Linear finite element solution
        mesh = genesis_mesh(os.path.join(data_path, "QuarterCylinderQuad4.g"))
        V = FEModel(mesh=mesh)
        el = Element(type="CPE4")
        mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
        V.assign_properties(
            element_block="ElementBlock1", element_type=el, material=mat
        )
        V.dirichlet_bc("Nodeset-200", X)
        V.dirichlet_bc("Nodeset-201", Y)

        step = V.static_step()
        # Pressure on inside face
        step.pressure("SURFACE-1", 1.0)
        step.run()
        filename = "VolumeLocking_%d.png" % i
        if plot:
            plot2d(
                model=V,
                deformed=1,
                color="b",
                linestyle="-.",
                label=r"Linear, $\nu=%g$" % Nu,
                ax=ax,
                # show=1,
                filename=filename,
                xlim=(-0.2, 5),
                ylim=(-0.2, 5),
            )


if __name__ == "__main__":
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_volume_locking(data_path, plot=True)
