import logging
from numpy import *
from numpy.linalg import det, inv

from ..x.constants import *
from ..x.utilities import *
from .element_base import element_base


class isop_base(element_base):
    """Base class for isoparametric stress-displacement elements"""

    num_gauss = None

    @classmethod
    def num_integration(klass):
        return klass.num_gauss

    @staticmethod
    def gauss_rule_info(point=None):
        raise NotImplementedError

    def bmatrix(self, *args):
        raise NotImplementedError

    def shapefun_der(self, *args):
        raise NotImplementedError

    @property
    def numdof(self):
        return sum([count_digits(nfs) for nfs in self.signature])

    @property
    def numdofpernod(self):
        return count_digits(self.signature[0])

    def pmatrix(self, N):
        n = count_digits(self.signature[0])
        P = zeros((n, self.nodes * n))
        for i in range(self.dimensions):
            P[i, i :: self.dimensions] = N
        return P

    @classmethod
    def interpolate_to_centroid(cls, data, index=None):
        """Inverse distance weighted average of integration point data at the
        element centroid"""
        if index is not None:
            ntens = cls.ndir + cls.nshr
            m = len(cls.variables()) * ntens
            data = row_stack(
                [
                    data[(m * p) + index * ntens : (m * p) + (index + 1) * ntens]
                    for p in range(cls.num_gauss)
                ]
            )
        return cls.average(cls.cp, data)

    @classmethod
    def project_to_nodes(cls, data, v):
        """Inverse distance weighted average of integration point data at each
        element node"""
        nx = len(v)
        a = zeros((cls.nodes, nx))
        for i in range(cls.nodes):
            a[i, :] = cls.average(cls.xp[i], data, v)
        return a

    @classmethod
    def average(cls, point, data, v=None):
        """Inverse distance weighted average of integration point data at point"""

        if data.ndim == 1:
            # SCALAR DATA
            assert len(data) == cls.num_gauss

        elif len(data.shape) == 2:
            # VECTOR OR TENSOR DATA
            assert data.shape[0] == cls.num_gauss

        else:
            raise TypeError("Unknown data type")

        if cls.num_gauss == 1:
            weights = [1.0]
        else:
            dist = lambda a, b: max(sqrt(dot(a - b, a - b)), 1e-6)
            weights = zeros(cls.num_gauss)
            for p in range(cls.num_gauss):
                xi, _ = cls.gauss_rule_info(p)
                weights[p] = 1.0 / dist(point, xi)

        if data.ndim == 1:
            # SCALAR DATA
            return average(data, weights=weights)

        elif len(data.shape) == 2:
            # VECTOR OR TENSOR DATA
            return average(data, axis=0, weights=weights)
