import numpy as np

from felab.util.lang import is_listlike


# ---------------------------------------------------------------------------- #
# --- GAUSS INTEGRATION RULES ------------------------------------------------ #
# ---------------------------------------------------------------------------- #
def line_gauss_rule_info(rule=2, p=None):
    """Gauss integration rules

    Parameters
    ----------
    rule : int
        The integration rule
        1 - single point integration
        2 - 2 point integration
        3 - 3 point integration
        4 - 4 point integration

    Returns
    -------
    qcoord : ndarray
        qcoord[i] is the ith Gauss point
    w : ndarray
        w[i] is the ith Gauss weight

    """
    RT3, RT5 = np.sqrt(3.0), np.sqrt(5.0)
    w14 = w44 = 0.5 - np.sqrt(5.0 / 6.0) / 6.0
    w24 = w34 = 0.5 + np.sqrt(5.0 / 6.0) / 6.0
    s34 = np.sqrt((3 - 2 * np.sqrt(6.0 / 5.0)) / 7.0)
    s24 = -s34
    s44 = np.sqrt((3 + 2 * np.sqrt(6.0 / 5.0)) / 7.0)
    s14 = -s44
    gw = {
        1: (np.zeros(1), np.array([2.0])),
        2: (np.array([-1.0, 1.0]) / RT3, np.ones(2)),
        3: (np.array([-RT3, 0.0, RT3]) / RT5, np.array([5.0, 8.0, 5.0]) / 9.0),
        4: (np.array([s14, s24, s34, s44]), np.array([w14, w24, w34, w44])),
    }
    g, w = gw[rule]
    if p is None:
        return g, w
    return g[p], w[p]


def tri_gauss_rule_info(rule, point):
    if rule == 1:
        gp, w = np.ones(3) / 3, 1
    elif rule == 3:
        gp, w = np.ones(3) / 6, 1.0 / 3.0
        gp[point] = 2.0 / 3.0
    elif rule == -3:
        gp, w = np.ones(3) / 2, 1.0 / 3.0
        gp[point] = 0
    elif rule == 6:
        g1 = (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2.0 / 5.0))) / 18.0
        g2 = (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2.0 / 5.0))) / 18.0
        if point < 3:
            gp = np.array([g1, g1, g1])
            w = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
            gp[point] = 1 - 2 * g1
        else:
            gp = np.array([g2, g2, g2])
            w = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
            gp[point] = 1 - 2 * g2
    elif rule == 7:
        g1 = (6 - np.sqrt(15)) / 21
        g2 = (6 + np.sqrt(15)) / 21
        if point < 3:
            gp = np.array([g1, g1, g1])
            w = (155 - np.sqrt(15)) / 1200
            gp[point] = 1 - 2 * g1
        elif point > 2 and point < 6:
            gp = np.array([g2, g2, g2])
            w = (155 + np.sqrt(15)) / 1200
            gp[point] = 1 - 2 * g2
        elif point == 6:
            gp, w = np.ones(3) / 3.0, 9.0 / 40.0
    else:
        raise ValueError("Incorrect rule {0}".format(rule))
    return gp, w


def quad_gauss_rule_info(rule=2, point=None):
    """Gauss integration rules

    Parameters
    ----------
    rule : int
        The integration rule

    Returns
    -------
    qcoord : ndarray
        qcoord[i] is the ith Gauss point
    w : ndarray
        w[i] is the ith Gauss weight

    """
    if is_listlike(rule):
        p1, p2 = rule[:2]
    else:
        p1 = p2 = rule

    if point is None:
        xi, w1 = line_gauss_rule_info(p1)
        eta, w2 = line_gauss_rule_info(p2)
        return xi, w1, eta, w2

    j = int(np.floor((point - 1) / p1)) + 1
    i = point - p1 * (j - 1)
    xi, w1 = line_gauss_rule_info(p1, i - 1)
    eta, w2 = line_gauss_rule_info(p2, j - 1)

    return np.array([xi, eta]), w1 * w2
