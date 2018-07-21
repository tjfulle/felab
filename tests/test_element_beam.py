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
    Ke = El.response(zeros((2,2)),zeros((2,2)),[0,0],1.,1,1,svars,[],[],[],
                     STATIC, False, STIFF_ONLY, DIRECT)
