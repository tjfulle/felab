import numpy as np
from felab import *
from felab.io.plot import plot2d

mu = 1.0
Nu = 0.499
E = 2.0 * mu * (1.0 + Nu)


def demo_linear(data_path, ax=None):
    mesh = abaqus_mesh(os.path.join(data_path, "ThickCylinder_Linear.inp"))
    V = FEModel(jobid="VolumeLocking.Linear", mesh=mesh)
    el = Element(type="CPE4R")
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(element_block="ALL", element_type=el, material=mat, t=1)
    V.dirichlet_bc("SymYZ", X)
    V.dirichlet_bc("SymXZ", Y)
    step = V.static_step()
    # Pressure on inside face
    step.pressure("SurfID", 1.0)
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


def demo_quadratic(data_path, ax=None):
    mesh = abaqus_mesh(os.path.join(data_path, "ThickCylinder_Quadratic.inp"))
    V = FEModel(jobid="VolumeLocking.Quadratic", mesh=mesh)
    el = Element(type="CPE8B")
    mat = Material(name="Material-1", elastic=dict(E=E, Nu=Nu))
    V.assign_properties(element_block="ALL", element_type=el, material=mat, t=1)
    V.dirichlet_bc("SymYZ", X)
    V.dirichlet_bc("SymXZ", Y)
    # Pressure on inside face
    step = V.static_step()
    step.pressure("SurfID", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        return ax
    return None


def demo_analytic(data_path, plot=False):
    mesh = Mesh(filename=os.path.join(data_path, "ThickCylinder_Linear.inp"))
    ix = np.where(mesh.coord[:, 1] <= 1e-12)
    a = mesh.coord[ix][:, 0].min()
    b = mesh.coord[ix][:, 0].max()
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
            xy=u + mesh.coord, ax=None, color="orange", weight=6, label="Analytic"
        )
        return ax
    return None


def runall(plot=False):
    import os
    import matplotlib.pyplot as plt
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")

    ax = demo_analytic(data_path, plot=plot)
    ax = demo_linear(data_path, ax)
    demo_quadratic(data_path, )

    plt.legend()
    plt.show()


if __name__ == "__main__":
    runall(plot=True)
