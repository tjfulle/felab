from numpy import *
from conf import *
from felab import *

def test_element_beam_0():
    """Beam-Column stiffness test"""
    coord = array([[0, 0], [3, 4]], dtype=float)
    class mat: E = 100
    svars = zeros((2, 3))
    A, Izz = 125, 250
    El = B2D2(1, [0, 1], coord, mat, A=A, Izz=Izz)
    rhs = zeros(6)
    A = zeros((6, 6))
    lflags = [None, None, STIFF_ONLY, None, None, None]
    El.response(rhs, A, [], None, zeros((2,2)), zeros((2,2)), None, None, [0,0],
                1., 1, 1, None, None, None, lflags, None, None, None)
