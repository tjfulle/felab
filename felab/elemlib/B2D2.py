import numpy as np

from felab.constants import SCALAR, STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY
from .element_base import element_base


class B2D2(element_base):
    """Base class for 2 node Euler-Bernouli beam elements

    Parameters
    ----------
    label : int
        Element label
    elenod : list of int
        Internal node IDs of nodes making up this element
    elecoord : ndarray
        Coordinates of nodes
    elemat : object
        A felab.mat.Material instance
    elefab : dict
        Requires area 'A' and 'Izz'

    """

    nodes = 2
    dimensions = 2
    signature = [(1, 1, 0, 0, 0, 1, 0), (1, 1, 0, 0, 0, 1, 0)]
    elefab = {"A": None, "Izz": None}

    @classmethod
    def variables(cls):
        return (("P", SCALAR), ("S", SCALAR))

    def eval(
        self,
        rhs,
        A,
        svars,
        energy,
        u,
        du,
        v,
        a,
        time,
        dtime,
        kstep,
        kframe,
        dltyp,
        dlmag,
        predef,
        lflags,
        ddlmag,
        mdload,
        pnewdt,
    ):

        if lflags[2] not in (STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY):
            raise NotImplementedError

        if kframe == 0:
            return

        if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
            rhs[:] = 0.0
            if lflags[2] == RHS_ONLY:
                return

        # COMPUTE ELEMENT NORMAL
        v = self.xc[1] - self.xc[0]
        h = np.sqrt(np.dot(v, v))
        n = v / h

        # TRANSFORMATION MATRIX
        Te = np.eye(6)
        Te[0:2, 0:2] = Te[3:5, 3:5] = [[n[0], n[1]], [-n[1], n[0]]]

        # COLUMN STIFFNESS
        EA, EI = self.material.E * self.A, self.material.E * self.Izz
        K1 = (
            EA
            / h
            * np.array(
                [
                    [1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                ]
            )
        )
        # BEAM STIFFNESS
        K2 = (
            2.0
            * EI
            / h ** 3
            * np.array(
                [
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0],
                    [0.0, 6.0, 3.0 * h, 0.0, -6.0, 3.0 * h],
                    [0.0, 3.0 * h, 2.0 * h * h, 0.0, -3.0 * h, h * h],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, -6.0, -3.0 * h, 0.0, 6.0, -3.0 * h],
                    [0.0, 3.0 * h, h * h, 0.0, -3.0 * h, 2.0 * h * h],
                ]
            )
        )

        A[:] = np.dot(np.dot(Te.T, K1 + K2), Te)
