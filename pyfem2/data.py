from numpy import *
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
        self.temp_frame = None

    def __len__(self):
        return len(self.frames)

    @property
    def last_converged_frame(self):
        for frame in reversed(self.frames):
            if frame.converged:
                return frame

    def Frame(self, dtime, copy=1):
        frame = Frame(self.time, dtime)
        self.time += dtime
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = deepcopy(frame_n.field_outputs)
        frame.number = len(self.frames)
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
        self.converged = False

    def adjust_dt(self, dtime):
        self.increment = dtime
        self.value = self.start + dtime

    def SymmetricTensorField(self, name, position, labels,
                             ndir, nshr, *args, **kwargs):
        field = SymmetricTensorField(name, position, labels,
                                     ndir, nshr, *args, **kwargs)
        if position in (INTEGRATION_POINT, ELEMENT_CENTROID):
            name = (args[-1], name)
        self.field_outputs[name] = field
        return field

    def VectorField(self, name, position, labels, ncomp, *args, **kwargs):
        field = VectorField(name, position, labels, ncomp, *args, **kwargs)
        if position in (INTEGRATION_POINT, ELEMENT_CENTROID):
            name = (args[-1], name)
        self.field_outputs[name] = field
        return field

    def ScalarField(self, name, position, labels, *args, **kwargs):
        field = ScalarField(name, position, labels, *args, **kwargs)
        if position in (INTEGRATION_POINT, ELEMENT_CENTROID):
            name = (args[-1], name)
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

class FieldOutputs(OrderedDict):
    def __getitem__(self, key):
        try:
            return super(FieldOutputs, self).__getitem__(key)
        except KeyError as E:
            pass
        keys = []
        for k in self.keys():
            if is_listlike(k) and k[1] == key:
                # element block property
                keys.append(k)
        if not keys:
            raise E
        a = self[keys[0]]
        for k in keys[1:]:
            a = row_stack((a, self[k]))
        return a

class FieldOutput(object):

    def __init__(self, name, position, labels, type, components, shape, eleblk,
                 **kwargs):
        self.name = name
        self.position = position
        self.labels = labels
        self.type = type
        self.components = components
        self.key = 'displ' if name == 'U' else name
        self.eleblk = eleblk
        self._values = None
        if components is not None:
            self.keys = [self.key+x for x in components]
        else:
            self.keys = [self.key]
        self.data = zeros(shape)

        a = kwargs.get('data')
        if a is not None:
            self.add_data(a)

    def __getitem__(self, i):
        return self.data[i]

    def add_data(self, data, ix=None):
        if not is_listlike(data):
            data = ones_like(self.data) * data
        else:
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
            if self.position == ELEMENT_CENTROID:
                return self.data
            if self.position != INTEGRATION_POINT:
                raise ValueError('Cannot project data to centroid')
            # Interpolate Gauss point data to element center
            return array([InterpolateToCentroid(x) for x in self.data])

        raise ValueError('Unknown position')

    @property
    def values(self):
        if self._values is not None:
            return self._values
        self._values = []
        for (i, label) in enumerate(self.labels):
            fv = FieldValue(self.position, label, self.type,
                            self.components, self.data[i])
            self._values.append(fv)
        return self._values

class SymmetricTensorField(FieldOutput):
    def __init__(self, name, position, labels, ndir, nshr, *args, **kwargs):
        eleblk = None
        components = ('xx', 'yy', 'zz')[:ndir] + ('xy', 'yz', 'xz')[:nshr]
        if position == INTEGRATION_POINT:
            ngauss, eleblk = args
        if ndir is not None and nshr is not None:
            ntens = ndir + nshr
        else:
            ntens = 0
        num = len(labels)
        if position == INTEGRATION_POINT:
            if ngauss:
                shape = (num, ngauss, ntens)
            else:
                shape = (num, ntens)
        else:
            shape = (num, ntens)
        super(SymmetricTensorField, self).__init__(
            name, position, labels, SYMTENSOR, components, shape, eleblk, **kwargs)

class VectorField(FieldOutput):
    def __init__(self, name, position, labels, nc, *args, **kwargs):
        eleblk = None
        num = len(labels)
        components = ('x', 'y', 'z')[:nc]
        if position == INTEGRATION_POINT:
            ngauss, eleblk = args
            shape = (num, ngauss, nc)
        else:
            shape = (num, nc)
        super(VectorField, self).__init__(
            name, position, labels, VECTOR, components, shape, eleblk, **kwargs)

class ScalarField(FieldOutput):
    def __init__(self, name, position, labels, *args, **kwargs):
        eleblk = None
        num = len(labels)
        components = None
        if position == INTEGRATION_POINT:
            ngauss, eleblk = args
            shape = (num, ngauss)
        else:
            shape = (num, 1)
        super(ScalarField, self).__init__(
            name, position, labels, SCALAR, components, shape, eleblk, **kwargs)

class FieldValue:
    def __init__(self, position, label, type, components, data):
        self.position = position
        self.label = label
        self.data = data
        self.type = type
        self.components = components
        self._mag = None
        self._prin = None

    @property
    def magnitude(self):
        if self._mag is None:
            if self.type in (SCALAR, VECTOR):
                self._mag = sqrt(dot(self.data, self.data))
            else:
                w = array([1.,1.,1.,2.,2.,2.])
                self._mag = sqrt(sum(self.data * self.data * w))
        return self._mag

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
    def prinvals(self):
        if self.type != SYMTENSOR:
            return None
        if self._prin is None:
            sx, sy, sz, txy, tyz, txz = self._tensor_components()
            s3, s2, s1 = sorted(eigvalsh([[sx,txy,txz],[txy,sy,tyz],[txz,tyz,sz]]))
            self._prin = array([s3, s2, s1])
        return self._prin

    @property
    def max_principal(self):
        if self.type != SYMTENSOR:
            return None
        return self.prinvals[-1]

    @property
    def min_principal(self):
        if self.type != SYMTENSOR:
            return None
        p = self.prinvals
        if len(self.components) == 6:
            return p[0]
        else:
            return p[1]
