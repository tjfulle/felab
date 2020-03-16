from numpy import *
from numpy.linalg import eigvalsh
from collections import OrderedDict

from ..x.utilities import *

__all__ = ["FieldOutputs", "ScalarField", "VectorField", "SymmetricTensorField"]


class FieldOutputs(OrderedDict):
    def __getitem__(self, key):
        try:
            return super(FieldOutputs, self).__getitem__(key)
        except KeyError as E:
            pass
        keys = []
        for k in self.keys():
            if is_listlike(k) and k[1] == key:
                # ELEMENT BLOCK PROPERTY
                keys.append(k)
        if not keys:
            raise KeyError(key)
        a = self[keys[0]]
        for k in keys[1:]:
            a = row_stack((a, self[k]))
        return a


class FieldOutput(object):
    def __init__(
        self,
        name,
        position,
        labels,
        type,
        components,
        shape,
        eleblk,
        elements=None,
        data=None,
    ):
        self.name = name
        self.position = position
        self.labels = labels
        self.type = type
        self.components = components
        self.key = "displ" if name == "U" else name
        self.eleblk = eleblk
        self._values = None
        if components is not None:
            self.keys = [self.key + x for x in components]
        else:
            self.keys = [self.key]
        self.shape = shape
        self.data = zeros(self.shape)

        if data is not None:
            idata = asarray(data)
            if len(idata) == shape[-1]:
                idata = zeros(shape)
                idata[:] = data
            self.add_data(idata)

        self._elements = elements
        reqels = position in (ELEMENT_CENTROID, INTEGRATION_POINT)
        if reqels:
            if elements is None:
                raise TypeError("Expected elements to be passed")
        else:
            if elements is not None:
                raise TypeError("Field position not element based")

    def __getitem__(self, i):
        return self.data[i]

    def add_data(self, data, ix=None):
        if self.data is None:
            self.data = zeros(self.shape)
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
                raise ValueError("Cannot project data to centroid")
            # INTERPOLATE GAUSS POINT DATA TO ELEMENT CENTER
            # FIXME: is the index i below right?  it used to be e
            return array(
                [
                    self._elements[i].interpolate_to_centroid(x)
                    for (i, x) in enumerate(self.data)
                ]
            )

        raise ValueError("Unknown position")

    @property
    def values(self):
        if self._values is not None:
            return self._values
        self._values = []
        for (i, label) in enumerate(self.labels):
            fv = FieldValue(
                self.position, label, self.type, self.components, self.data[i]
            )
            self._values.append(fv)
        return self._values


class SymmetricTensorField(FieldOutput):
    def __init__(
        self,
        name,
        position,
        labels,
        ndir,
        nshr,
        eleblk=None,
        ngauss=None,
        elements=None,
        data=None,
    ):

        if position == INTEGRATION_POINT and not ngauss:
            raise TypeError("Expected ngauss")

        components = ("xx", "yy", "zz")[:ndir] + ("xy", "yz", "xz")[:nshr]

        num = len(labels)
        ntens = ndir + nshr
        if ngauss:
            shape = (num, ngauss, ntens)
        else:
            shape = (num, ntens)

        super(SymmetricTensorField, self).__init__(
            name,
            position,
            labels,
            SYMTENSOR,
            components,
            shape,
            eleblk,
            elements=elements,
            data=data,
        )


class TensorField(FieldOutput):
    def __init__(
        self,
        name,
        position,
        labels,
        ndir,
        nshr,
        eleblk=None,
        ngauss=None,
        elements=None,
        data=None,
    ):

        if position == INTEGRATION_POINT and not ngauss:
            raise TypeError("Expected ngauss")

        components = array([["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]])
        if ndir < 3:
            components[2, 2] = ""
            if ndir < 2:
                components[1, 1] = ""
        if nshr < 3:
            components[0, 2] = components[2, 0] = ""
            if nshr < 2:
                components[1, 2] = components[2, 1] = ""
        components = [x for x in components.flatten() if x.split()]

        num = len(labels)
        ntens = ndir + 2 * nshr
        if ngauss:
            shape = (num, ngauss, ntens)
        else:
            shape = (num, ntens)

        super(TensorField, self).__init__(
            name,
            position,
            labels,
            TENSOR,
            components,
            shape,
            eleblk,
            elements=elements,
            data=data,
        )


class VectorField(FieldOutput):
    def __init__(
        self,
        name,
        position,
        labels,
        nc,
        eleblk=None,
        ngauss=None,
        elements=None,
        data=None,
    ):

        if position == INTEGRATION_POINT and not ngauss:
            raise TypeError("Expected ngauss")

        components = ("x", "y", "z")[:nc]

        num = len(labels)
        if ngauss:
            shape = (num, ngauss, nc)
        else:
            shape = (num, nc)

        super(VectorField, self).__init__(
            name,
            position,
            labels,
            VECTOR,
            components,
            shape,
            eleblk,
            elements=elements,
            data=data,
        )


class ScalarField(FieldOutput):
    def __init__(
        self, name, position, labels, eleblk=None, ngauss=None, elements=None, data=None
    ):

        if position == INTEGRATION_POINT and not ngauss:
            raise TypeError("Expected ngauss")

        components = None

        num = len(labels)
        if ngauss:
            shape = (num, ngauss)
        else:
            shape = (num, 1)

        super(ScalarField, self).__init__(
            name,
            position,
            labels,
            SCALAR,
            components,
            shape,
            eleblk,
            elements=elements,
            data=data,
        )


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
                w = array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])
                self._mag = sqrt(sum(self.data * self.data * w))
        return self._mag

    def _tensor_components(self):
        if len(self.components) == 3:
            # PLANE STRESS
            sx, sy, txy = self.data
            sz, tyz, txz = 0, 0, 0
        elif len(self.components) == 4:
            # PLANE STRAIN
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
            s3, s2, s1 = sorted(
                eigvalsh([[sx, txy, txz], [txy, sy, tyz], [txz, tyz, sz]])
            )
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
