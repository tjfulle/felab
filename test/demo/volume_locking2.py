import numpy as np
from felab import *
from felab.io.plot import plot2d

mu = 1.0
Nu = 0.499
E = 2.0 * mu * (1.0 + Nu)


def demo_linear(data_path, ax=None):
    mesh = genesis_mesh(os.path.join(data_path, "QuarterCylinderQuad4.g"))
    V = FEModel(jobid="VolumeLocking.Linear", mesh=mesh)

    el = Element(type="CPE4")
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
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


def demo_bbar(data_path, ax=None):
    mesh = genesis_mesh(os.path.join(data_path, "QuarterCylinderQuad4.g"))
    V = FEModel(jobid="VolumeLocking.BBar", mesh=mesh)
    el = Element(type="CPE4B")
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
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


def demo_reduced_integration(data_path, ax=None):
    mesh = genesis_mesh(os.path.join(data_path, "QuarterCylinderQuad4.g"))
    V = FEModel(jobid="VolumeLocking.Reduced", mesh=mesh)
    el = Element(type="CPE4R")
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
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


def demo_selectively_reduced_integration(data_path, ax=None):
    mesh = genesis_mesh(os.path.join(data_path, "QuarterCylinderQuad4.g"))
    V = FEModel(jobid="VolumeLocking.SelReduced", mesh=mesh)
    el = Element(type="CPE4RS")
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
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


def demo_quadratic(data_path, ax=None):
    mesh = genesis_mesh(os.path.join(data_path, "QuarterCylinderQuad8.g"))
    V = FEModel(jobid="VolumeLocking.Quadratic", mesh=mesh)
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    el = Element(type="CPE8B")
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
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


def demo_analytic(data_path, plot=False):
    mesh = Mesh(filename=os.path.join(data_path, "QuarterCylinderQuad4.g"))
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
    import os
    import matplotlib.pyplot as plt
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")

    ax = demo_analytic(data_path, plot=True)
    ax = demo_reduced_integration(data_path, ax)
    ax = demo_selectively_reduced_integration(data_path, ax)
    ax = demo_linear(data_path, ax)
    demo_quadratic(data_path)

    plt.legend()
    plt.show()


if __name__ == "__main__":
    runall(plot=True)
