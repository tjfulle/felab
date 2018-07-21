from numpy import *
from conf import *
from felab import *

def test_element_elastic_link_0():
    """ElasticLink2 stiffness test"""
    class mat: E = 1
    svars = zeros((2, 3))
    El = L1D2(1, [0, 1], [0, 1], mat, A=1)
    Ke = El.response(zeros(2), zeros(2), [0,0], 1., 1, 1, svars, [], [], STATIC,
                     [], False, STIFF_ONLY, DIRECT)
    assert allclose([[1,-1],[-1,1]], Ke)

def test_element_elastic_link_1():
    """ElasticLink2 stiffness test"""
    class mat: E=1000
    svars = zeros((2, 3))
    El = L2D2(1, [0, 1], [[0,0], [30,40]], mat, A=5)
    Ke = El.response(zeros((2,2)),zeros((2,2)),[0,0],1.,1,1,svars,[],[],STATIC,
                     [], False, STIFF_ONLY, DIRECT)
    assert allclose([[ 36.,  48., -36., -48.],
                     [ 48.,  64., -48., -64.],
                     [-36., -48.,  36.,  48.],
                     [-48., -64.,  48.,  64.]], Ke)

def test_element_elastic_link_2():
    """ElasticLink2 stiffness test"""
    class mat: E = 343
    svars = zeros((2, 3))
    El = L3D2(1, [0, 1], [[0,0,0],[2,3,6]], mat, A=10)
    Ke = El.response(zeros((2,3)),zeros((2,3)),[0,0],1.,1,1,svars,[],[],STATIC,
                     [], False, STIFF_ONLY, DIRECT)
    assert allclose([[  40.,   60.,  120.,  -40.,  -60., -120.],
                     [  60.,   90.,  180.,  -60.,  -90., -180.],
                     [ 120.,  180.,  360., -120., -180., -360.],
                     [ -40.,  -60., -120.,   40.,   60.,  120.],
                     [ -60.,  -90., -180.,   60.,   90.,  180.],
                     [-120., -180., -360.,  120.,  180.,  360.]], Ke)

