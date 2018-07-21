import logging
from numpy import *
from ..utilities import UserInputError, is_listlike
from .elas import elas
from .linear_elastic import linear_elastic
from .neohooke import neo_hooke
from .thermal import thermally_conductive

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
         mat.elastic(E=10e6, Nu=.29)

    - Create a material with isotropic heat conduction :math:`\kappa=12`:

      .. code:: python

         mat = Material('Material-1')
         mat.isotropic_thermal_conductivity(12)

    """
    def __init__(self, name, **kwds):

        self.name = name
        self._model = None

        # YOUNG'S MODULUS AND POISSON'S RATIO
        self.E, self.Nu = None, None

        # THERMAL CONDUCTIVITY
        self.k_iso = None

        # DENSITY
        self.density = None

        # CHECK KEYWORDS
        for (kwd, v) in kwds.items():
            k = kwd.lower()
            if k in ('density', 'rho'):
                self.density = v

            elif k == 'elastic':
                try:
                    self.elastic(**v)
                except TypeError:
                    raise UserInputError('ELASTIC PROPERTIES MUST BE A '
                                         'DICT, NOT {0}'.format(type(v)))

            elif 'neo' in k and 'hooke' in k:
                try:
                    self.neo_hooke(**v)
                except TypeError:
                    raise UserInputError('NEOHOOKE PROPERTIES MUST BE A '
                                         'DICT, NOT {0}'.format(type(v)))

            elif 'thermal_conductivity' in k:
                self.thermal_conductivity(v)

            else:
                logging.warn('SETTING UNKNOWN MATERIAL PROPERTY {0!r}'.format(kwd))
                setattr(self, kwd, v)

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, m):
        self._model = m
        if not hasattr(self._model, 'requires'):
            self._model.requires = None

    def Density(self, rho):
        """Assign mass density

        Parameters
        ----------
        rho : float
            Mass density

        """
        assert rho > 0
        self.density = rho

    def elastic(self, **kwds):
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

             self.elastic(E=10, Nu=.29)

        - Bulk modulus :math:`K=10` and shear modulus :math:`G=6`:

          .. code:: python

             self.elastic(K=10, G=6)

        """
        self.set_elastic_properties(**kwds)
        self.model = linear_elastic(self.Lame, self.G)

    def neo_hooke(self, **kwds):
        """Assign Neo Hooke properties

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

             self.neo_hooke(E=10, Nu=.29)

        - Bulk modulus :math:`K=10` and shear modulus :math:`G=6`:

          .. code:: python

             self.neo_hooke(K=10, G=6)

        """
        self.set_elastic_properties(**kwds)
        self.model = neo_hooke(self.E, self.Nu)

    def set_elastic_properties(self, **kwds):
        if 'rho' in [s.lower() for s in kwds.keys()]:
            if len(kwds) != 3:
                raise UserInputError('EXACTLY 2 ELASTIC CONSTANTS REQUIRED')
            for (k,v) in kwds.items():
                if k.lower() == 'rho':
                    self.Density(v)
        elif len(kwds) != 2:
            raise UserInputError('EXACTLY 2 ELASTIC CONSTANTS REQUIRED')
        props = elas(**kwds)
        # ADD PROPERTIES TO SELF
        self.__dict__.update(props)

    def thermal_conductivity(self, k):
        """Assign the coefficient of thermal conductivity

        Parameters
        ----------
        k : ndarray or float
            Coefficient of thermal conductivity

        """
        self.model = thermally_conductive(k)
    isotropic_thermal_conductivity = thermal_conductivity

    def response(self, stress, statev, strain, dstrain, time, dtime,
                 temp, dtemp, predef, dpred, ndir, nshr, ntens,
                 coords, F0, F, noel, kstep, kinc):
        return self.model.response(
            stress, statev, strain, dstrain, time, dtime, temp, dtemp,
            predef, dpred, ndir, nshr, ntens, coords, F0, F, noel, kstep, kinc)
