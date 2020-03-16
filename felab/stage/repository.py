from felab.stage.stage import load_stage
from felab.stage.diffusive_ht import diffusive_ht_stage
from felab.stage.static import static_stage
from felab.stage.dynamic import dynamic_stage

__all__ = ["repository"]


class repository(object):
    def __init__(self, model):
        self.model = model
        self._keys = []
        self._values = []

    def __setitem__(self, key, value):
        self._keys.append(key)
        self._values.append(value)

    def __getitem__(self, key):
        return self._values[self._keys.index(key)]
        try:
            i = self._keys.index(key)
        except ValueError:
            raise KeyError(key)
        return self._values[i]

    def __len__(self):
        return len(self._keys)

    def __iter__(self):
        return iter(self._keys)

    def __contains__(self, key):
        return key in self._keys

    def keys(self):
        return self._keys

    def values(self):
        return self._values

    def items(self):
        for (i, key) in enumerate(self._keys):
            yield (key, self._values[i])

    def InitialStage(self, name):
        period = 0.0
        stage = load_stage(self.model, len(self), name, None, period)
        self[name] = stage
        return self.last

    def create_static_stage(self, name, period=1.0, copy=True, **kwds):
        last = self._values[-1].increments[-1]
        if not last.converged:
            raise RuntimeError("PREVIOUS STAGE HAS UNCONVERGED INCREMENTS")
        stage = static_stage(self.model, len(self), name, self.last, period, **kwds)
        if copy:
            stage.copy_from(self.last)
        stage.increments[0].converged = True
        self[name] = stage
        return self.last

    def create_dynamic_stage(self, name, period=1.0, copy=True, **kwds):
        last = self._values[-1].increments[-1]
        if not last.converged:
            raise RuntimeError("PREVIOUS STAGE HAS UNCONVERGED INCREMENTS")
        stage = dynamic_stage(self.model, len(self), name, self.last, period, **kwds)
        if copy:
            stage.copy_from(self.last)
        stage.increments[0].converged = True
        self[name] = stage
        return self.last

    def create_heat_transfer_stage(self, name, period, copy=True):
        last = self._values[-1].increments[-1]
        if not last.converged:
            raise RuntimeError("PREVIOUS STAGE HAS UNCONVERGED INCREMENTS")
        stage = diffusive_ht_stage(self.model, len(self), name, self.last, period)
        if copy:
            stage.copy_from(self.last)
        stage.increments[0].converged = True
        self[name] = stage
        return self.last

    @property
    def last(self):
        return self._values[-1]

    @property
    def first(self):
        return self._values[0]
