from numpy import *
from copy import deepcopy

from ..utilities import *
from ..constants import *
from ..data_wharehouse import *

__all__ = ['StepRepository1']

# --------------------------------------------------------------------------- #
# --- THE FOLLOWING IS FOR OUTPUT ONLY -------------------------------------- #
# --------------------------------------------------------------------------- #
class StepRepository1(object):
    def __init__(self):
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
    def keys(self):
        return self._keys
    def values(self):
        return self._values
    def items(self):
        for (i, key) in enumerate(self._keys):
            yield (key, self._values[i])
    def Step(self, name=None, copy=1):
        if name is None:
            j = 1
            while 1:
                name = 'Step-{0}'.format(len(self._keys)+j)
                if name not in self._keys:
                    break
                j += 1
        if not self._values:
            step = Step(name, 0.)
        else:
            last = self._values[-1].framces[-1]
            step = Step(name, last.value)
            if copy:
                step.frames[0].field_outputs = deepcopy(last.field_outputs)
        self[name] = step
        return self._values[-1]
    @property
    def last(self):
        return self._values[-1]
    @property
    def first(self):
        return self._values[0]

class Step(object):
    def __init__(self, name, start):
        self.written = 0
        self.name = name
        self.time = start
        self.frames = []
        self.Frame(0.)

    def __len__(self):
        return len(self.frames)

    def Frame(self, dtime, copy=1):
        frame = Frame(self.time, dtime)
        self.time += dtime
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = deepcopy(frame_n.field_outputs)
        frame.number = len(self.frames)
        self.frames.append(frame)
        return frame

class Frame(object):
    def __init__(self, start, dtime):
        self.start = start
        self.increment = dtime
        self.value = start + dtime
        self.field_outputs = FieldOutputs()
        self.converged = False

    def adjust_dt(self, dtime):
        self.increment = dtime
        self.value = self.start + dtime

    def SymmetricTensorField(self, name, position, labels, ndir, nshr,
                             eleblk=None, ngauss=None,  elements=None, data=None):
        field = SymmetricTensorField(name, position, labels, ndir, nshr,
                                     eleblk=eleblk, ngauss=ngauss,
                                     elements=elements, data=data)

        if field.eleblk is not None:
            name = (field.eleblk, name)
        self.field_outputs[name] = field

        return field

    def VectorField(self, name, position, labels, ncomp, eleblk=None,
                    ngauss=None, elements=None, data=None):

        field = VectorField(name, position, labels, ncomp, eleblk=eleblk,
                            ngauss=ngauss, elements=elements, data=data)

        if field.eleblk is not None:
            name = (field.eleblk, name)
        self.field_outputs[name] = field

        return field

    def ScalarField(self, name, position, labels, eleblk=None, ngauss=None,
                    elements=None, data=None):

        field = ScalarField(name, position, labels, eleblk=eleblk,
                            ngauss=ngauss, elements=elements, data=data)

        if field.eleblk is not None:
            name = (field.eleblk, name)
        self.field_outputs[name] = field

        return field

    def add_data(self, **kwds):
        for (key, value) in kwds.items():
            d = {}
            if isinstance(value, tuple):
                if len(value) == 2:
                    d['ix'] = value[1]
                    value = value[0]
                else:
                    raise ValueError('Unknown add_data option for {0}'.format(key))
            self.field_outputs[key].add_data(value, **d)
