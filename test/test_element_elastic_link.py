from numpy import *
from conf import *
from felab import *

lflags = [None, None, STIFF_ONLY, None, None, None]


def test_element_elastic_link_0():
    """ElasticLink2 stiffness test"""

    class mat:
        E = 1

    svars = zeros((2, 3))
    El = L1D2(1, [0, 1], [0, 1], mat, A=1)
    rhs = zeros(2)
    A = zeros((2, 2))
    El.response(
        rhs,
        A,
        svars,
        None,
        zeros(2),
        zeros(2),
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
    assert allclose([[1, -1], [-1, 1]], A)


def test_element_elastic_link_1():
    """ElasticLink2 stiffness test"""

    class mat:
        E = 1000

    svars = zeros((2, 3))
    El = L2D2(1, [0, 1], [[0, 0], [30, 40]], mat, A=5)
    rhs = zeros(4)
    A = zeros((4, 4))
    El.response(
        rhs,
        A,
        svars,
        None,
        zeros((2, 2)),
        zeros((2, 2)),
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
    assert allclose(
        [
            [36.0, 48.0, -36.0, -48.0],
            [48.0, 64.0, -48.0, -64.0],
            [-36.0, -48.0, 36.0, 48.0],
            [-48.0, -64.0, 48.0, 64.0],
        ],
        A,
    )


def test_element_elastic_link_2():
    """ElasticLink2 stiffness test"""

    class mat:
        E = 343

    svars = zeros((2, 3))
    El = L3D2(1, [0, 1], [[0, 0, 0], [2, 3, 6]], mat, A=10)
    rhs = zeros(6)
    A = zeros((6, 6))
    El.response(
        rhs,
        A,
        svars,
        None,
        zeros((2, 3)),
        zeros((2, 3)),
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
    assert allclose(
        [
            [40.0, 60.0, 120.0, -40.0, -60.0, -120.0],
            [60.0, 90.0, 180.0, -60.0, -90.0, -180.0],
            [120.0, 180.0, 360.0, -120.0, -180.0, -360.0],
            [-40.0, -60.0, -120.0, 40.0, 60.0, 120.0],
            [-60.0, -90.0, -180.0, 60.0, 90.0, 180.0],
            [-120.0, -180.0, -360.0, 120.0, 180.0, 360.0],
        ],
        A,
    )
