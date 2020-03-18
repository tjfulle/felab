import numpy as np

from felab.error import UserInputError
from felab.util.lang import is_listlike
from felab.material.material import MaterialModel


class ThermallyConductive(MaterialModel):
    name = "Thermally Conductive"

    def __init__(self, *, k):
        if not is_listlike(k):
            self.k_iso = k
        else:
            k = np.asarray(k)
            self.k_aniso = np.eye(3)
            if len(k) == 3:
                np.fill_diagonal(self.k_aniso, k)
            elif k.size == 9:
                self.k_aniso[:] = k.reshape(3, 3)
            else:
                raise UserInputError("K MUST BE A 3 VECTOR OR 3X3 ARRAY")
            self.k_iso = np.trace(self.k_aniso) / 3.0

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
        predef,
        dpred,
        ndir,
        nshr,
        ntens,
        coords,
        F0,
        F,
        noel,
        kstep,
        kframe,
    ):
        return self.isotropic_thermal_conductivity(coords.shape[0])

    def isotropic_thermal_conductivity(self, ndim):
        return self.k_iso * np.eye(ndim)
