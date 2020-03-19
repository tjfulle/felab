from felab.error import UserInputError
from felab.material.elas import elas


__all__ = ["Material", "MaterialModel"]


class Material:
    """The base material object.

    Parameters
    ----------
    name : str
        The name of the material
    kwds : dict
        Keywords to be assigned as material attributes.

    """  # noqa: W605

    def __init__(self, *, name, model, density=None):
        self.name = name
        self.model = model
        self.density = density
        self.E = getattr(model, "E", None)
        self.G = getattr(model, "G", None)

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, arg):
        if arg is not None:
            assert arg > 0
        self._density = arg

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, arg):
        assert isinstance(arg, MaterialModel)
        self._model = arg

    def eval(self, *args):
        return self.model.eval(*args)


class MaterialModel(object):
    requires = None

    def compute_elastic_properties(self, **kwds):
        if "rho" in [s.lower() for s in kwds.keys()]:
            if len(kwds) != 3:
                raise UserInputError("EXACTLY 2 ELASTIC CONSTANTS REQUIRED")
            for (k, v) in kwds.items():
                if k.lower() == "rho":
                    self.Density(v)
        elif len(kwds) != 2:
            raise UserInputError("EXACTLY 2 ELASTIC CONSTANTS REQUIRED")
        props = elas(**kwds)
        return props

    def eval(self, *args):
        """Evaluate the material model"""
        raise NotImplementedError
