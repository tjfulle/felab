import numpy as np
from felab import *
from felab.elemlib import CPE8B, CPE4R, CPE4B, CPE4RS, CPE4
from felab.io.plot import plot2d

mu = 1.0
Nu = 0.499
E = 2.0 * mu * (1.0 + Nu)


def demo_linear(ax=None):
    mesh = genesis_mesh("./data/QuarterCylinderQuad4.g")
    V = FEModel(jobid="VolumeLocking.Linear", mesh=mesh)

    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE4, material=mat, t=1
    )

    V.dirichlet_bc("Nodeset-200", X)
    V.dirichlet_bc("Nodeset-201", Y)
    # Pressure on inside face
    step = V.static_step()
    step.pressure("SURFACE-1", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        ax = plot2d(
            model=V,
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


def demo_bbar(ax=None):
    mesh = genesis_mesh("./data/QuarterCylinderQuad4.g")
    V = FEModel(jobid="VolumeLocking.BBar", mesh=mesh)
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE4B, material=mat, t=1
    )
    V.dirichlet_bc("Nodeset-200", X)
    V.dirichlet_bc("Nodeset-201", Y)
    # Pressure on inside face
    step = V.static_step()
    step.pressure("SURFACE-1", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        ax = plot2d(model=V, deformed=1, ax=ax, color="r", label="bbar")
        return ax
    return None


def demo_reduced_integration(ax=None):
    mesh = genesis_mesh("./data/QuarterCylinderQuad4.g")
    V = FEModel(jobid="VolumeLocking.Reduced", mesh=mesh)
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE4R, material=mat, t=1
    )
    V.dirichlet_bc("Nodeset-200", X)
    V.dirichlet_bc("Nodeset-201", Y)
    # Pressure on inside face
    step = V.static_step()
    step.pressure("SURFACE-1", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        ax = plot2d(model=V, deformed=1, ax=ax, color="k", label="reduced")
        return ax
    return None


def demo_selectively_reduced_integration(ax=None):
    mesh = genesis_mesh("./data/QuarterCylinderQuad4.g")
    V = FEModel(jobid="VolumeLocking.SelReduced", mesh=mesh)
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE4RS, material=mat, t=1
    )
    V.dirichlet_bc("Nodeset-200", X)
    V.dirichlet_bc("Nodeset-201", Y)
    # Pressure on inside face
    step = V.static_step()
    step.pressure("SURFACE-1", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        ax = plot2d(
            model=V, deformed=1, ax=ax, color="g", label="Sel. reduced integration"
        )
        return ax
    return None


def demo_quadratic(ax=None):
    mesh = genesis_mesh("./data/QuarterCylinderQuad8.g")
    V = FEModel(jobid="VolumeLocking.Quadratic", mesh=mesh)
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=CPE8B, material=mat, t=1
    )
    V.dirichlet_bc("Nodeset-200", X)
    V.dirichlet_bc("Nodeset-201", Y)
    # Pressure on inside face
    step = V.static_step(solver=NEWTON)
    step.pressure("SURFACE-1", 1.0)
    # V.surface_load("SURFACE-300", [0.195090322, 0.98078528])
    # V.surface_load("SURFACE-301", [0.555570233, 0.831469612])
    # V.surface_load("SURFACE-302", [0.831469612, 0.555570233])
    # V.surface_load("SURFACE-303", [0.98078528, 0.195090322])
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
        ax = plot2d(
            mesh=mesh,
            xy=u + mesh.coord,
            ax=None,
            color="orange",
            weight=8,
            label="Analytic",
        )
        return ax
    return None


def runall(plot=True):
    import matplotlib.pyplot as plt

    ax = demo_analytic(plot=True)
    ax = demo_reduced_integration(ax)
    ax = demo_selectively_reduced_integration(ax)
    ax = demo_linear(ax)
    demo_quadratic()

    plt.legend()
    plt.show()


if __name__ == "__main__":
    runall(plot=True)
