import numpy as np
from felab.fe_model import fe_model
from felab.mesh import Mesh
from felab.elemlib import CPE4
from felab.constants import X, Y


def demo_volume_locking(plot=False):
    mu = 1.0
    for (i, Nu) in enumerate((0.0, 0.2, 0.45, 0.499)):
        E = 2.0 * mu * (1.0 + Nu)

        # Analytic solution
        mesh = Mesh(filename="./data/QuarterCylinderQuad4.g")
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
            ax = mesh.Plot2D(
                xy=u + mesh.coord,
                color="orange",
                label=r"Analytic, $\nu=%g$" % Nu,
                weight=8,
            )

        # Linear finite element solution
        V = fe_model()
        V.genesis_mesh("./data/QuarterCylinderQuad4.g")
        mat = V.material("Material-1")
        mat.elastic(E=E, Nu=Nu)
        V.assign_properties("ElementBlock1", CPE4, mat)
        V.assign_prescribed_bc("Nodeset-200", X)
        V.assign_prescribed_bc("Nodeset-201", Y)

        step = V.static_step()
        # Pressure on inside face
        step.assign_pressure("SURFACE-1", 1.0)
        step.run()
        filename = "VolumeLocking_%d.png" % i
        if plot:
            V.Plot2D(
                deformed=1,
                color="b",
                linestyle="-.",
                label=r"Linear, $\nu=%g$" % Nu,
                ax=ax,
                #show=1,
                filename=filename,
                xlim=(-0.2, 5),
                ylim=(-0.2, 5),
            )


if __name__ == "__main__":
    demo_volume_locking(plot=True)
