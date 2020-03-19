import numpy as np
from numpy.linalg import inv

from felab.material.material import MaterialModel


class LinearElastic(MaterialModel):
    """Linear elastic material

    Parameters
    ----------
    density : float
        The mass density
    kwds : dict
        name=value pairs of elastic constants

    Notes
    -----
    Any two elastic constants can be passed to this method from which all
    other elastic constants are calculated (isotropy of the tangent
    elastic stiffness assumed).

    Examples
    --------
    - Young's modulus :math:`E=10` and Poisson's ratio :math:`\\nu=.29`:

      .. code:: python

         LinearElastic(E=10, Nu=.29)

    - Bulk modulus :math:`K=10` and shear modulus :math:`G=6`:

      .. code:: python

         LinearElastic(K=10, G=6)

    """

    name = "Linear Elastic"

    def __init__(self, **kwds):
        props = self.compute_elastic_properties(**kwds)
        self.Lambda = props["Lame"]
        self.Mu = props["G"]
        self.__dict__.update(**props)

    def eval(
        self,
        stress,
        statev,
        strain,
        dstrain,
        time,
        dtime,
        temp,
        dtemp,
        ndir,
        nshr,
        ntens,
    ):

        C11 = self.Lambda + 2 * self.Mu
        C12 = self.Lambda
        C44 = self.Mu
        D = np.array(
            [
                [C11, C12, C12, 0, 0, 0],
                [C12, C11, C12, 0, 0, 0],
                [C12, C12, C11, 0, 0, 0],
                [0, 0, 0, C44, 0, 0],
                [0, 0, 0, 0, C44, 0],
                [0, 0, 0, 0, 0, C44],
            ]
        )

        if nshr == 1:
            # MODIFY THE STIFFNESS FOR 2D ACCORDING TO:
            # 1) PLANE STRAIN: REMOVE ROWS AND COLUMNS OF THE STIFFNESS
            #    CORRESPONDING TO THE PLANE OF ZERO STRAIN
            # 2) PLANE STRESS: INVERT THE STIFFNESS AND REMOVE THE ROWS
            #    AND COLUMNS OF THE COMPLIANCE CORRESPONDING THE PLANE OF
            #    ZERO STRESS.
            if ndir == 2:
                # PLANE STRESS
                # INVERT THE STIFFNESS TO GET THE COMPLIANCE
                idx = [[[0], [1], [3]], [0, 1, 3]]
                D = inv(inv(D)[idx])
            elif ndir == 3:
                # PLANE STRAIN
                idx = [[[0], [1], [2], [3]], [0, 1, 2, 3]]
                D = D[idx]

        stress += np.dot(D, dstrain)
        return stress, statev, D
