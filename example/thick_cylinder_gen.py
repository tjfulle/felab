import numpy as np
from felab.fe_model import fe_model
from felab.elemlib import CPE4, CPE4R, CPE4RS, CPE8B
from felab.constants import X, Y
from felab.mesh import Mesh

mu = 1.0
Nu = 0.499
E = 2.0 * mu * (1.0 + Nu)


def demo_linear(ax=None):
    V = fe_model(jobid="VolumeLocking.Linear")
    V.genesis_mesh("./data/QuarterCylinderQuad4.g")
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties("ElementBlock1", CPE4, mat, t=1)
    V.assign_prescribed_bc("Nodeset-200", X)
    V.assign_prescribed_bc("Nodeset-201", Y)
    # Pressure on inside face
    step = V.static_step()
    step.assign_pressure("SURFACE-1", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        ax = V.Plot2D(
            deformed=1,
            color="b",
            linestyle="-.",
            ax=ax,
            label="Linear",
            xlim=(-0.2, 5),
            ylim=(-0.2, 5),
        )
        return ax
    return None


def demo_reduced_integration(ax=None):
    V = fe_model(jobid="VolumeLocking.Reduced")
    V.genesis_mesh("./data/QuarterCylinderQuad4.g")
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties("ElementBlock1", CPE4R, mat, t=1, hourglass_control=True)
    V.assign_prescribed_bc("Nodeset-200", X)
    V.assign_prescribed_bc("Nodeset-201", Y)

    step = V.static_step()
    # Pressure on inside face
    step.assign_pressure("SURFACE-1", 1.0)
    step.run()

    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, ax=ax, color="b", label="Reduced integration")
        return ax
    return None


def demo_selectively_reduced_integration(ax=None):
    V = fe_model(jobid="VolumeLocking.SelReduced")
    V.genesis_mesh("./data/QuarterCylinderQuad4.g")
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties("ElementBlock1", CPE4RS, mat, t=1)
    V.assign_prescribed_bc("Nodeset-200", X)
    V.assign_prescribed_bc("Nodeset-201", Y)

    step = V.static_step()
    # Pressure on inside face
    step.assign_pressure("SURFACE-1", 1.0)
    step.run()

    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, ax=ax, color="b", label="Sel reduced integration")
        return ax
    return None


def demo_quadratic(ax=None):
    V = fe_model(jobid="VolumeLocking.Quadratic")
    V.genesis_mesh("./data/QuarterCylinderQuad8.g")
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties("ElementBlock1", CPE8B, mat, t=1)
    V.assign_prescribed_bc("Nodeset-200", X)
    V.assign_prescribed_bc("Nodeset-201", Y)

    step = V.static_step()
    # Pressure on inside face
    # V.assign_pressure('SURFACE-1', 1.)
    step.assign_surface_load("SURFACE-300", [0.195090322, 0.98078528])
    step.assign_surface_load("SURFACE-301", [0.555570233, 0.831469612])
    step.assign_surface_load("SURFACE-302", [0.831469612, 0.555570233])
    step.assign_surface_load("SURFACE-303", [0.98078528, 0.195090322])
    step.run()

    V.write_results()


def demo_analytic(plot=False):
    mesh = Mesh(filename="./data/QuarterCylinderQuad4.g")
    a = mesh.coord[0, 1]
    b = mesh.coord[-1, 0]
    p = 1.0
    u = np.zeros_like(mesh.coord)
    for (i, x) in enumerate(mesh.coord):
        r = np.sqrt(x[0] ** 2 + x[1] ** 2)
        term1 = (1.0 + Nu) * a ** 2 * b ** 2 * p / (E * (b ** 2 - a ** 2))
        term2 = 1.0 / r + (1.0 - 2.0 * Nu) * r / b ** 2
        ur = term1 * term2
        u[i, :] = ur * x[:] / r
    mesh.put_nodal_solution("VolumeLocking.Analytic", u)
    if plot:
        ax = mesh.Plot2D(
            xy=u + mesh.coord, ax=None, color="orange", weight=8, label="Analytic"
        )
        return ax
    return None


def runall(plot=False):
    import matplotlib.pyplot as plt

    ax = demo_analytic(plot=plot)
    # ax = demo_reduced_integration(ax)
    ax = demo_selectively_reduced_integration(ax)
    ax = demo_linear(ax)
    demo_quadratic()

    plt.legend()
    plt.show()


if __name__ == "__main__":
    runall(plot=True)
