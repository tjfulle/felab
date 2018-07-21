from numpy import *

from .utilities import *
from .constants import *
from .step import load_step
from .heat_transfer_step import heat_transfer_step
from .static_step import static_step
from .dynamic_step import dynamic_step

__all__ = ['step_repository']

class step_repository(object):
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

    def InitialStep(self, name):
        start, period = 0., 0.
        step = load_step(self.model, len(self), name, None, period)
        self[name] = step
        return self.last

    def create_static_step(self, name, period=1., **kwds):
        last = self._values[-1].frames[-1]
        if not last.converged:
            raise RuntimeError('PREVIOUS STEP HAS UNCONVERGED FRAMES')
        step = static_step(self.model, len(self), name, self.last, period, **kwds)
        if copy:
            step.copy_from(self.last)
        step.frames[0].converged = True
        self[name] = step
        return self.last

    def create_dynamic_step(self, name, period=1., **kwds):
        last = self._values[-1].frames[-1]
        if not last.converged:
            raise RuntimeError('PREVIOUS STEP HAS UNCONVERGED FRAMES')
        step = dynamic_step(self.model, len(self), name, self.last, period, **kwds)
        if copy:
            step.copy_from(self.last)
        step.frames[0].converged = True
        self[name] = step
        return self.last

    def create_heat_transfer_step(self, name, period):
        last = self._values[-1].frames[-1]
        if not last.converged:
            raise RuntimeError('PREVIOUS STEP HAS UNCONVERGED FRAMES')
        step = heat_transfer_step(self.model, len(self), name, self.last, period)
        if copy:
            step.copy_from(self.last)
        step.frames[0].converged = True
        self[name] = step
        return self.last

    @property
    def last(self):
        return self._values[-1]

    @property
    def first(self):
        return self._values[0]
