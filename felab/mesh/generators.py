import os
import numpy as np

from felab.mesh import Mesh
from felab.error import UserInputError


__all__ = [
    "unit_square_mesh",
    "rectilinear_mesh2d",
    "vtk_mesh",
    "abaqus_mesh",
    "genesis_mesh",
]


def rectilinear_mesh2d(
    nx=1, ny=1, lx=1, ly=1, shiftx=None, shifty=None, method=None, order=1
):
    """Generates a rectilinear 2D finite element mesh.

    Parameters
    ----------
    nx, ny : int
        number of elements in x and y directions, respectively
    lx, ly : float
        length of mesh bounding box in x and y directions, respectively

    Returns
    -------
    Mesh object

    """
    if nx < 1:
        raise UserInputError("Requres at least 1 element in x")
    if ny < 1:
        raise UserInputError("Requres at least 1 element in y")
    if lx < 0.0:
        raise UserInputError("Requres positive length in x")
    if ly < 0.0:
        raise UserInputError("Requres positive length in y")

    if order == 1:
        if method in (2, "felippa"):
            p, t = gen_quad4_pt_2(nx, ny, lx, ly, sx=shiftx, sy=shifty)
        elif method in (None, "default", 1):
            p, t = gen_quad4_pt_1(nx, ny, lx, ly, sx=shiftx, sy=shifty)
        else:
            raise UserInputError("Wrong 2D meshing method for 1st order mesh")
    elif order == 2:
        if method not in (None, 2, "felippa"):
            raise UserInputError("Wrong 2D meshing method for 2nd order mesh")
        p, t = gen_quad8_pt_2(nx, ny, lx, ly, sx=shiftx, sy=shifty)
    else:
        raise UserInputError("Unsupported 2D meshing order")

    nodtab, eletab = gen_node_and_elem_tables(p, t)
    return Mesh(nodtab=nodtab, eletab=eletab)


def unit_square_mesh(nx=1, ny=1, shiftx=None, shifty=None, method=None, order=1):
    """
    Generates a rectilinear 2D finite element mesh.

    Parameters
    ----------
    shape : tuple
        (nx, ny) where nx is the number elements in :math:`x` and ny
        number of element in :math:`y`.

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls the ``Mesh.rectilinear_mesh_2d`` class method and
    stores the returned mesh as the ``fe_model.mesh`` attribute.

    """
    return rectilinear_mesh2d(
        nx=nx,
        ny=ny,
        lx=1,
        ly=1,
        shiftx=shiftx,
        shifty=shifty,
        method=method,
        order=order,
    )


def gen_node_and_elem_tables(p, t):
    nodtab = [[n + 1] + list(x) for (n, x) in enumerate(p)]
    eletab = [[e + 1] + list(c + 1) for (e, c) in enumerate(t)]
    return nodtab, eletab


def gen_quad4_pt_1(nx, ny, lx, ly, sx=None, sy=None):

    sx, sy = sx or 0.0, sy or 0.0
    xpoints = np.linspace(0, lx, nx + 1) + sx
    ypoints = np.linspace(0, ly, ny + 1) + sy
    p = np.array([(x, y) for y in ypoints for x in xpoints])

    # Connectivity
    k = 0
    numele = nx * ny
    t = np.zeros((numele, 4), dtype=int)
    for elem_num in range(numele):
        ii = elem_num + k
        elem_nodes = [ii, ii + 1, ii + nx + 2, ii + nx + 1]
        t[elem_num, :] = elem_nodes
        if (elem_num + 1) % (nx) == 0:
            k += 1

    return p, t


def gen_quad4_pt_2(nx, ny, lx, ly, sx=None, sy=None):
    """Generate a rectilinear mesh with quad4 elements"""

    # Bounding box
    sx, sy = sx or 0.0, sy or 0.0
    xc = np.array(
        [
            [0.0 + sx, 0.0 + sy],
            [lx + sx, 0.0 + sy],
            [lx + sx, ly + sy],
            [0.0 + sx, ly + sy],
        ]
    )

    def shape(x, y):
        a = np.array(
            [
                (1.0 - x) * (1.0 - y),
                (1.0 + x) * (1.0 - y),
                (1.0 + x) * (1.0 + y),
                (1.0 - x) * (1.0 + y),
            ]
        )
        return a / 4.0

    # Nodal coordinates
    numnod = (nx + 1) * (ny + 1)
    p = np.zeros((numnod, 2))
    k = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            xi = 2.0 * i / nx - 1.0
            eta = 2.0 * j / ny - 1.0
            N = shape(xi, eta)
            p[k] = (np.dot(N, xc[:, 0]), np.dot(N, xc[:, 1]))
            k += 1

    # Element connectivity
    numele = nx * ny
    t = np.zeros((numele, 4), dtype=int)
    k = 0
    for i in range(nx):
        for j in range(ny):
            c1 = (ny + 1) * i + j
            c2 = (c1 + 1) + ny
            t[k] = (c1, c2, c2 + 1, c1 + 1)
            k += 1

    return p, t


def gen_quad8_pt_2(nx, ny, lx, ly, sx=None, sy=None):
    """Generate a rectilinear mesh with quad8 elements"""

    # Bounding box
    sx, sy = sx or 0.0, sy or 0.0
    xc = np.array([[sx, sy], [lx + sy, sy], [lx + sx, ly + sy], [sx, ly + sy]])

    def shape(x, y):
        a = np.array(
            [
                (1.0 - x) * (1.0 - y),
                (1.0 + x) * (1.0 - y),
                (1.0 + x) * (1.0 + y),
                (1.0 - x) * (1.0 + y),
            ]
        )
        return a / 4.0

    # Nodal coordinates
    numnod = (2 * nx + 1) * (2 * ny + 1) - nx * ny
    p = np.zeros((numnod, 2))
    k = 0
    for i in range(2 * nx + 1):
        for j in range(2 * ny + 1):
            if (i + 1) % 2 == 0 and (j + 1) % 2 == 0:
                continue
            xi = float(i - nx) / nx
            eta = float(j - ny) / ny
            N = shape(xi, eta)
            p[k] = (np.dot(N, xc[:, 0]), np.dot(N, xc[:, 1]))
            k += 1

    # Element connectivity
    numele = nx * ny
    t = np.zeros((numele, 8), dtype=int)
    k = 0
    for i in range(nx):
        for j in range(ny):
            c1 = (3 * ny + 2) * i + 2 * (j + 1) - 2
            c3 = c1 + 2 * ny - j + 1
            c2 = c3 + ny + j + 1
            t[k] = (c1, c2, c2 + 2, c1 + 2, c3, c2 + 1, c3 + 1, c1 + 1)
            k += 1

    return p, t


def genesis_mesh(filename):
    """
    Generates a finite element mesh from a Genesis file.

    Parameters
    ----------
    filename : str
        The path to a valid Genesis file

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls ``mesh.Mesh`` with the ``filename`` keyword and
    stores the returned mesh as the ``fe_model.mesh`` attribute.

    """
    if not os.path.isfile(filename):
        raise UserInputError("NO SUCH FILE {0!r}".format(filename))
    return Mesh(filename=filename)


def abaqus_mesh(filename):
    """
    Generates a finite element mesh from a Abaqus input file.

    Parameters
    ----------
    filename : str
        The path to a valid Genesis file

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls ``mesh.Mesh`` with the ``filename`` keyword and
    stores the returned mesh as the ``fe_model.mesh`` attribute.

    """
    if not os.path.isfile(filename):
        raise UserInputError("NO SUCH FILE {0!r}".format(filename))
    return Mesh(filename=filename)


def vtk_mesh(filename):
    """
    Generates a finite element mesh from a vtk .vtu file.

    Parameters
    ----------
    filename : str
        The path to a valid .vtu file

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls ``mesh.Mesh`` with the ``filename`` keyword and
    stores the returned mesh as the ``fe_model.mesh`` attribute.

    """
    if not os.path.isfile(filename):
        raise UserInputError("NO SUCH FILE {0!r}".format(filename))
    return Mesh(filename=filename)
