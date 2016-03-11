from numpy import *
import logging
from numpy.linalg import eigvalsh
from collections import OrderedDict
from copy import deepcopy

from .utilities import *
from .elemlib1 import InterpolateToCentroid

__all__ = ['StepRepository']

class StepRepository(object):
    def __init__(self):
        self._keys = []
        self._values = []
    def __setitem__(self, key, value):
        self._keys.append(key)
        self._values.append(value)
    def __getitem__(self, key):
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
        self.temp_frame = None

    def __len__(self):
        return len(self.frames)

    def Frame(self, dtime, copy=1):
        frame = Frame(self.time, dtime)
        self.time += dtime
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = deepcopy(frame_n.field_outputs)
        self.frames.append(frame)
        return frame

    def TempFrame(self):
        self.temp_frame = Frame(self.time, 0.)
        self.temp_frame.field_outputs = deepcopy(self.frames[-1].field_outputs)

class Frame(object):
    def __init__(self, start, dtime):
        self.start = start
        self.increment = dtime
        self.value = start + dtime
        self.field_outputs = FieldOutputs()

    def adjust_dt(self, dtime):
        self.increment = dtime
        self.value = self.start + dtime

    def SymmetricTensorField(self, name, position, labels, ndir, nshr, *args):
        field = SymmetricTensorField(name, position, labels, ndir, nshr, *args)
        if position in (INTEGRATION_POINT, ELEMENT):
            name = (args[1], name)
        self.field_outputs[name] = field

    def VectorField(self, name, position, labels, ncomp, *args):
        field = VectorField(name, position, labels, ncomp, *args)
        if position in (INTEGRATION_POINT, ELEMENT):
            name = (args[1], name)
        self.field_outputs[name] = field

    def ScalarField(self, name, position, labels, *args):
        field = ScalarField(name, position, labels, *args)
        if position in (INTEGRATION_POINT, ELEMENT):
            name = (args[1], name)
        self.field_outputs[name] = field

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

class FieldOutputs(OrderedDict):
    pass

class FieldOutput(object):

    def __init__(self, name, position, labels, type, components=None):
        self.name = name
        self.position = position
        self.labels = labels
        self.type = type
        self.comopnents = components
        self.key = 'displ' if name == 'U' else name
        if components is not None:
            self.keys = [self.key+x for x in components]
        else:
            self.keys = [self.key]

    def __getitem__(self, i):
        return self.data[i]

    def add_data(self, data, ix=None):
        data = asarray(data)
        if ix is not None:
            self.data[ix] = data
        else:
            assert data.size == self.data.size
            self.data[:] = reshape(data, self.data.shape)

    def get_data(self, position=None):

        if position is None:
            return self.data

        if position == ELEMENT_CENTROID:

            if self.position != INTEGRATION_POINT:
                raise ValueError('Cannot project data to centroid')

            # Interpolate Gauss point data to element center
            return array([InterpolateToCentroid(x) for x in self.data])

        raise ValueError('Unknown position')

    @property
    def values(self):
        raise NotImplementedError
        if hasattr(self, '_values'):
            return self._values
        self._values = []
        for (i, label) in enumerate(self.labels):
            fv = FieldValue(self.position, label, self.type,
                            self.components, self.data[i])
            self._values.append(fv)
        return self._values

class SymmetricTensorField(FieldOutput):
    def __init__(self, name, position, labels, ndir, nshr, *args):
        c = ('xx', 'yy', 'zz')[:ndir] + ('xy', 'yz', 'xz')[:nshr]
        super(SymmetricTensorField, self).__init__(name, position, labels,
                                                   SYMTENSOR, c)
        if self.position == INTEGRATION_POINT:
            ngauss, eleblk = args
            self.eleblk = eleblk
        if ndir is not None and nshr is not None:
            ntens = ndir + nshr
        else:
            ntens = 0
        num = len(self.labels)
        if position == INTEGRATION_POINT:
            if ngauss:
                shape = (num, ngauss, ntens)
            else:
                shape = (num, ntens)
        else:
            shape = (num, ntens)
        self.data = zeros(shape)

class VectorField(FieldOutput):
    def __init__(self, name, position, labels, nc, *args):
        c = ('x', 'y', 'z')[:nc]
        super(VectorField, self).__init__(name, position, labels, VECTOR, c)
        num = len(self.labels)
        if self.position == INTEGRATION_POINT:
            ngauss, self.eleblk = args
            shape = (num, ngauss, nc)
        else:
            shape = (num, nc)
        self.data = zeros(shape)

class ScalarField(FieldOutput):
    def __init__(self, name, position, labels, *args):
        super(ScalarField, self).__init__(name, position, labels, SCALAR)
        num = len(self.labels)
        if self.position == INTEGRATION_POINT:
            ngauss, self.eleblk = args
            shape = (num, ngauss)
        else:
            shape = (num, 1)
        self.data = zeros(shape)

class FieldValue:
    def __init__(self, position, label, type, components, data):
        self.position = position
        self.label = label
        self.data = data
        self.type = type
        self.components = components

    @property
    def magnitude(self):
        if self.type in (SCALAR, VECTOR):
            return sqrt(dot(self.data, self.data))


    def _tensor_components(self):
        if len(self.components) == 3:
            # Plane stress
            sx, sy, txy = self.data
            sz, tyz, txz = 0, 0, 0
        elif len(self.components) == 4:
            # Plane strain
            sx, sy, sz, txy = self.data
            tyz, txz = 0, 0
        else:
            sx, sy, sz, txy, tyz, txz = self.data
        return sx, sy, sz, txy, tyz, txz

    @property
    def max_principal(self):
        if self.type != SYMTENSOR:
            return None
        sx, sy, sz, txy, tyz, txz = self._tensor_components()
        s3, s2, s1 = sorted(eigvalsh([[sx,txy,txz],[txy,sy,tyz],[txz,tyz,sz]]))
        return s1

    @property
    def min_principal(self):
        if self.type != SYMTENSOR:
            return None
        sx, sy, sz, txy, tyz, txz = self._tensor_components()
        s3, s2, s1 = sorted(eigvalsh([[sx,txy,txz],[txy,sy,tyz],[txz,tyz,sz]]))
        if len(self.components) == 6:
            return s3
        else:
            return s2
