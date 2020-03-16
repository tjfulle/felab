from felab.first_fe_program import UniformBar


def test_model_uniform_bar_0():
    xa, xb = 0.0, 1.0
    A, E, p, q = 1, 1, 1, 0
    u, e, f = UniformBar(xa, xb, A, E, p, q, numele=10)
    assert abs(u[-1] - 1.0) < 1e-12


def test_model_uniform_bar_1():
    xa, xb = 0.0, 1.0
    A, E, p, q = 1, 1, 1, 1
    u, e, f = UniformBar(xa, xb, A, E, p, q, numele=1000)
    assert abs(u[-1] - 1.5) < 5e-4
