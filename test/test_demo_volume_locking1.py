from felab import *


def Runall(plot=False):
    mu = 1.0
    for (i, Nu) in enumerate((0.0, 0.2, 0.45, 0.499)):
        E = 2.0 * mu * (1.0 + Nu)

        filename = "../../felab-docs/_images/VolumeLocking_%d.png" % i

        # Analytic solution
        mesh = Mesh(filename="./data/QuarterCylinderQuad4.g")
        p = 1.0
        a, b = mesh.coord[0, 1], mesh.coord[-1, 0]
        u = zeros_like(mesh.coord)
        for (i, x) in enumerate(mesh.coord):
            r = sqrt(x[0] ** 2 + x[1] ** 2)
            term1 = (1.0 + Nu) * a ** 2 * b ** 2 * p / (E * (b ** 2 - a ** 2))
            term2 = 1.0 / r + (1.0 - 2.0 * Nu) * r / b ** 2
            ur = term1 * term2
            u[i, :] = ur * x[:] / r
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
        mat = V.create_material("Material-1")
        mat.elastic(E=E, Nu=Nu)
        V.assign_properties("ElementBlock1", CPE4, mat)
        V.assign_prescribed_bc("Nodeset-200", X)
        V.assign_prescribed_bc("Nodeset-201", Y)

        stage = V.create_static_stage()
        # Pressure on inside face
        stage.assign_pressure("SURFACE-1", 1.0)
        stage.run()
        if plot:
            V.Plot2D(
                deformed=1,
                color="b",
                linestyle="-.",
                label=r"Linear, $\nu=%g$" % Nu,
                ax=ax,
                show=1,  # filename=filename,
                xlim=(-0.2, 5),
                ylim=(-0.2, 5),
            )


def test(plot=False):
    if plot:
        import matplotlib.pyplot as plt

        plt.clf()
        plt.cla()
    Runall(plot=plot)


if __name__ == "__main__":
    test(plot=True)
