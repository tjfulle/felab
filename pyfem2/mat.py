import logging
from numpy import *
from numpy.linalg import inv
from .elas import elas
from .utilities import UserInputError, is_listlike

__all__ = ['Material']

class Material(object):
    """The base material object.

    Parameters
    ----------
    name : str
        The name of the material
    kwds : dict
        Keywords to be assigned as material attributes.

    Notes
    -----
    This is the main material class. Once instantiated, properties are
    assigned to the material by different methods.

    Examples
    --------
    - Create a linear elastic material with Young's modulus
      :math:`E=10\\times 10^6` and Poisson's ratio :math:`\\nu=.29`

      .. code:: python

         mat = Material('Material-1')
         mat.Elastic(E=10e6, Nu=.29)

    - Create a material with isotropic heat conduction :math:`\kappa=12`:

      .. code:: python

         mat = Material('Material-1')
         mat.IsotropicThermalConductivity(12)

    """
    def __init__(self, name, **kwds):
        self.name = name

        # Young's modulus and Poisson's ratio
        self.E, self.Nu = None, None

        # Thermal conductivity
        self.k_iso = None

        # Density
        self.density = None
        for (kwd, v) in kwds.items():
            k = kwd.lower()
            if k in ('density', 'rho'):
                self.density = v
            elif k == 'elastic':
                try:
                    self.Elastic(**v)
                except TypeError:
                    raise UserInputError('ELASTIC PROPERTIES MUST BE A '
                                         'DICT, NOT {0}'.format(type(v)))
            elif 'thermal_conductivity' in k:
                self.ThermalConductivity(v)
            else:
                logging.warn('SETTING UNKNOWN MATERIAL PROPERTY {0!r}'.format(kwd))
                setattr(self, kwd, v)

    def Density(self, rho):
        """Assign mass density

        Parameters
        ----------
        rho : float
            Mass density

        """
        assert rho > 0
        self.density = rho

    def Elastic(self, **kwds):
        """Assign elastic properties

        Parameters
        ----------
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

             self.Elastic(E=10, Nu=.29)

        - Bulk modulus :math:`K=10` and shear modulus :math:`G=6`:

          .. code:: python

             self.Elastic(K=10, G=6)

        """
        if 'rho' in [s.lower() for s in kwds.keys()]:
            if len(kwds) != 3:
                raise UserInputError('EXACTLY 2 ELASTIC CONSTANTS REQUIRED')
            for (k,v) in kwds.items():
                if k.lower() == 'rho':
                    self.Density(v)
        elif len(kwds) != 2:
            raise UserInputError('EXACTLY 2 ELASTIC CONSTANTS REQUIRED')
        props = elas(**kwds)
        self.E, self.Nu = props['E'], props['Nu']
        self.G, self.K = props['G'], props['K']
        self.Mu = self.G
        self.Lambda = props['Lame']

    def ThermalConductivity(self, k):
        """Assign the coefficient of thermal conductivity

        Parameters
        ----------
        k : ndarray or float
            Coefficient of thermal conductivity

        """
        if not is_listlike(k):
            self.k_iso = k
        else:
            k = asarray(k)
            self.k_aniso = eye(3)
            if len(k) == 3:
                fill_diagonal(self.k_aniso, k)
            elif k.size == 9:
                self.k_aniso[:] = k.reshape(3,3)
            else:
                raise UserInputError('K MUST BE A 3 VECTOR OR 3X3 ARRAY')
            self.k_iso = trace(self.k_aniso) / 3.
    IsotropicThermalConductivity = ThermalConductivity

    def isotropic_thermal_conductivity(self, ndim):
        """The isotropic thermal conductivity matrix"""
        return self.k_iso * eye(ndim)

    def response(self, stress, statev, strain, dstrain, time, dtime,
                 temp, dtemp, ndir, nshr):
        D = self.stiffness(ndir, nshr)
        stress = dot(D, strain+dstrain)
        return stress, statev, D

    def stiffness(self, ndir, nshr, disp=None):
        if self.E is None:
            raise UserInputError('ELASTIC MODULUS NOT SET')
        if self.Nu is None:
            raise UserInputError("POISSON'S RATIO NOT SET")
        D = self.isotropic_elastic_stiffness()

        if nshr == 1:
            # Modify the stiffness for 2D according to:
            # 1) Plane strain: Remove rows and columns of the stiffness
            #    corresponding to the plane of zero strain
            # 2) Plane stress: Invert the stiffness and remove the rows
            #    and columns of the compliance corresponding the plane of
            #    zero stress.
            if ndir == 2:
                # plane stress
                # Invert the stiffness to get the compliance
                idx = [[[0], [1], [3]], [0, 1, 3]]
                D = inv(inv(D)[idx])
            elif ndir == 3:
                # plane strain
                idx = [[[0], [1], [2], [3]], [0, 1, 2, 3]]
                D = D[idx]

        if disp is None:
            return D

        ntens = ndir + nshr
        if disp == 2:
            D1, D2 = zeros((ntens, ntens)), eye(ntens)
            D1[:ndir,:ndir] = D[0,1]
            return D1, D-D1

        raise NotImplementedError

    def isotropic_elastic_stiffness(self):
        """Tangent elastic stiffness"""
        c11 = self.Lambda + 2*self.Mu
        c12 = self.Lambda
        c44 = self.Mu
        C = array([[c11, c12, c12, 0,   0,   0  ],
                   [c12, c11, c12, 0,   0,   0  ],
                   [c12, c12, c11, 0,   0,   0  ],
                   [0,   0,   0,   c44, 0,   0  ],
                   [0,   0,   0,   0,   c44, 0  ],
                   [0,   0,   0,   0,   0,   c44]])
        return C
