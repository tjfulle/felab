from numpy import allclose, array
from felab.fe_model import fe_model
from felab.constants import ALL, Rr, Zr
from felab.elemlib import CAX4


def analytic(a, b, E, Nu, r, p):
    numer = a ** 2 * (1 + Nu) * (b ** 2 + r ** 2 * (1 - 2 * Nu))
    denom = E * (b ** 2 - a ** 2) * r
    ur = p * numer / denom
    return ur


def test_thick_pressurized_cylinder_quad4(plot=False):
    Kfac = 1.0
    nx, ny = 4, 1
    a, b, d, p = 4, 10, 2, 10
    E, Nu = 1000.0, 0.0
    V = fe_model()
    V.rectilinear_mesh(nx, ny, b - a, d, shiftx=a, method=2)
    mat = V.create_material("Material-1")
    mat.elastic(E=E, Nu=Nu)

    V.create_element_block("EALL", ALL)
    V.assign_properties("EALL", CAX4, mat, formulation=1)

    V.assign_prescribed_bc(ALL, Zr)

    stage = V.create_static_stage()
    pfor = Kfac * p * a * d
    stage.assign_concentrated_load((1, 2), Rr, pfor / 2.0)
    stage.run()

    if plot:
        ax = V.Plot2D(deformed=1, color="orange", weight=3)
        V.Plot2D(color="blue", weight=0.75, ax=ax, show=1)

    u = stage.increments[-1].field_outputs["U"]
    print("node model   analytic")
    for (i, urz) in enumerate(u.data):
        ur = urz[0]
        r = V.mesh.coord[i, 0]
        ur_a = analytic(a, b, E, Nu, r, p)
        print("{0:4d} {1:.4f} {2:.4f}".format(i + 1, ur, ur_a))

    # For the quad element, the results are not exact for 4 elements, thus, the
    # expected solution comparison (and not the analytic solution)
    expected = array(
        [
            0.054635477465204724,
            0.05463547746520472,
            0.04472949240053272,
            0.044729492400532714,
            0.04024936095840438,
            0.04024936095840437,
            0.038341661490036316,
            0.03834166149003633,
            0.03785419098608201,
            0.03785419098608201,
        ]
    )
    assert allclose(u.data[:, 0], expected)
    assert allclose(u.data[:, 1], 0)


if __name__ == "__main__":
    test_thick_pressurized_cylinder_quad4(plot=True)
