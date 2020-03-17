from numpy import zeros_like
from felab.fe_model import FEModel
from felab.elemlib import CPS4, CPS4I
from felab.constants import ALL, X, Y, ILO, IHI
from felab.io.plot import plot2d

mu = 10000.0
nu = 0.0
E = 2.0 * mu * (1.0 + nu)


def demo_plane_stress_beam(ratio=0.25, plot=False):
    V = FEModel()
    length = 10.0
    a = ratio * length
    P = 2.22 / length ** 3 * E * a ** 3
    q = P / a
    V.rectilinear_mesh(nx=10, ny=3, lx=length, ly=2 * a)
    mat = V.material("Material-1")
    mat.elastic(E=E, Nu=nu)
    V.element_block(name="ElementBlock1", elements=ALL)
    El = CPS4I
    El = CPS4
    V.assign_properties(
        element_block="ElementBlock1", element_type=El, material=mat, t=1
    )
    V.dirichlet_bc(IHI, (X, Y))

    step = V.static_step()
    step.surface_load(ILO, [0, -q])
    step.run()

    if plot:
        ax = plot2d(model=V, deformed=1, color="b")
    u = analytic(V.mesh.coord, q)
    xy = V.mesh.coord + u
    if plot:
        plot2d(mesh=V.mesh, xy=xy, color="r", ax=ax, show=1)


def analytic(coord, q):
    a = (coord[:, 1].max() - coord[:, 1].min()) / 2.0
    L = coord[:, 0].max() - coord[:, 0].min()
    b = 1.0
    P = 2 * q * a * b
    # d = -P * L ** 3 / (2 * E * a ** 3 * b)
    II = b * (2.0 * a) ** 3 / 12.0
    u = zeros_like(coord)
    for (i, x) in enumerate(coord):
        x1, x2 = x
        x2 -= a
        u[i, 0] = (
            P / (2.0 * E * II) * x1 ** 2 * x2
            + nu * P / (6 * E * II) * x2 ** 3
            - P / (6 * II * mu) * x2 ** 3
            - (P * L ** 2 / (2 * E * II) - P * a ** 2 / (2 * II * mu)) * x2
        )
        u[i, 1] = (
            -nu * P / (2 * E * II) * x1 * x2 ** 2
            - P / (6 * E * II) * x1 ** 3
            + P * L ** 2 / (2 * E * II) * x1
            - P * L ** 3 / (3 * E * II)
        )
    return u


if __name__ == "__main__":
    # ---- Change the aspect ratio to see shear locking
    ratio = 0.25
    demo_plane_stress_beam(ratio=ratio, plot=True)
