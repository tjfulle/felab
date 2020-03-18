import numpy as np
from felab.elemlib.B2D2 import B2D2
from felab.constants import STIFF_ONLY


def test_element_beam_0():
    """Beam-Column stiffness test"""
    coord = np.array([[0, 0], [3, 4]], dtype=float)

    class mat:
        E = 100

    A, Izz = 125, 250
    El = B2D2(1, [0, 1], coord, mat, A=A, Izz=Izz)
    rhs = np.zeros(6)
    A = np.zeros((6, 6))
    lflags = [None, None, STIFF_ONLY, None, None, None]
    El.eval(
        rhs,
        A,
        [],
        None,
        np.zeros((2, 2)),
        np.zeros((2, 2)),
        None,
        None,
        [0, 0],
        1.0,
        1,
        1,
        None,
        None,
        None,
        lflags,
        None,
        None,
        None,
    )
