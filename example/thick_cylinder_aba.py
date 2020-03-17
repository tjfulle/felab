import numpy as np
from felab.fe_model import fe_model
from felab.elemlib import CPE4R, CPE8B
from felab.constants import X, Y
from felab.mesh import Mesh
from felab.io.plot import plot2d

mu = 1.0
Nu = 0.499
E = 2.0 * mu * (1.0 + Nu)


def demo_linear(ax=None):
    V = fe_model(jobid="VolumeLocking.Linear")
    V.abaqus_mesh("./data/ThickCylinder_Linear.inp")
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties("ALL", CPE4R, mat, t=1)
    V.assign_prescribed_bc("SymYZ", X)
    V.assign_prescribed_bc("SymXZ", Y)
    step = V.static_step()
    # Pressure on inside face
    step.assign_pressure("SurfID", 1.0)
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


def demo_quadratic(ax=None):
    V = fe_model(jobid="VolumeLocking.Quadratic")
    V.abaqus_mesh("./data/ThickCylinder_Quadratic.inp")
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties("ALL", CPE8B, mat, t=1)
    V.assign_prescribed_bc("SymYZ", X)
    V.assign_prescribed_bc("SymXZ", Y)
    # Pressure on inside face
    step = V.static_step()
    step.assign_pressure("SurfID", 1.0)
    step.run()
    V.write_results()
    if ax is not None:
        return ax
    return None


def demo_analytic(plot=False):
    mesh = Mesh(filename="./data/ThickCylinder_Linear.inp")
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
    import matplotlib.pyplot as plt

    ax = demo_analytic(plot=plot)
    ax = demo_linear(ax)
    demo_quadratic()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    runall(plot=True)
