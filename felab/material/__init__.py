from felab.material.material import Material as _Material
from felab.material.neohooke import NeoHooke
from felab.material.linear_elastic import LinearElastic
from felab.material.thermal import ThermallyConductive


def Material(*, name=None, density=None, **kwds):
    if len(kwds) != 1:
        raise ValueError("One, and only one, model should be specified (for now)")
    for (key, val) in kwds.items():
        if key == "elastic":
            factory = LinearElastic
        elif key == "neohooke":
            factory = NeoHooke
        elif key in ("thermal_conductivity", "isotropic_thermal_conductivity"):
            if not isinstance(val, dict):
                val = {"k": val}
            factory = ThermallyConductive
        else:
            raise ValueError(f"Unknown model {key}")
        model = factory(**val)
        break
    return _Material(name=name, density=density, model=model)
