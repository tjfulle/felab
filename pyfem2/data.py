from numpy import *
from numpy.linalg import eigvalsh
import logging
from collections import OrderedDict

from .utilities import *

def component_labels(type, numdim, ndir, nshr):
    if type == SYMTENSOR:
        components = ['xx','yy','zz'][:ndir] + ['xy','yz','xz'][:nshr]
    elif type == VECTOR:
        components = ['x','y','z'][:numdim]
    else:
        components = None
    return components

class TimeStepRepository(object):

    def __init__(self, mesh, first=1):
        self.mesh = mesh
        self.time = 0.
        self.steps = []
        if first:
            self.TimeStep(0.)
    def __getitem__(self, i):
        return self.steps[i]
    def __iter__(self):
        return iter(self.steps)
    def __len__(self):
        return len(self.steps)

    def TimeStep(self, dtime, copy=0):
        ts = TimeStep(self.mesh, self.time, dtime)
        self.time += dtime
        if copy:
            tsn = self.steps[-1]
            for fo in tsn.field_outputs.values():
                fo1 = ts.FieldOutput(fo.name, fo.type, fo.position)
                data = tsn.field_outputs[fo.name].get_data()
                if fo1.position == NODE:
                    fo1.add_data(data)
                else:
                    for (i, dx) in enumerate(data):
                        fo1.add_data(dx, block=i)
        self.steps.append(ts)
        return ts

class TimeStep(object):
    def __init__(self, mesh, time, dtime, copy=None):
        self.mesh = mesh
        self.start = time
        self.increment = dtime
        self.value = time + dtime
        self.field_outputs = FieldOutputs()

    def FieldOutput(self, name, type, position, numcomp=None):
        fo = FieldOutput(self.mesh, name, type, position, numcomp=numcomp)
        self.field_outputs[name] = fo
        return fo

    def add_data(self, **kwargs):
        for (key, value) in kwargs.items():
            kwds = {}
            if isinstance(value, tuple):
                if len(value) == 2:
                    value = value[0]
                    kwds['ix'] = value[1]
                else:
                    raise ValueError('Unknown add_data option for {0}'.format(key))
            self.field_outputs[key].add_data(value, **kwds)

class FieldOutputs(OrderedDict):
    pass

class FieldOutput(object):
    def __init__(self, mesh, name, type, position, numcomp=None):
        self._values = None
        self.name = name
        self.key = 'displ' if name.upper() == 'U' else name
        self.type = type
        self.position = position
        self.mesh = mesh
        numdim, numele, numnod = mesh.numdim, mesh.numele, mesh.numnod

        if self.position not in (NODE, INTEGRATION_POINT, ELEMENT):
            raise ValueError('Unknown field position')

        if self.type not in (SCALAR, VECTOR, SYMTENSOR):
            raise ValueError('Unknown field type')

        if self.position == NODE:
            if self.type == SCALAR:
                self.components = None
                self.data = zeros((mesh.numnod,1))
            elif self.type == VECTOR:
                self.components = ('x','y','z')[:numdim]
                self.data = zeros((mesh.numnod, numdim))
            else:
                raise ValueError('Node tensors not yet supported')

        elif self.position == ELEMENT:
            if self.type == SCALAR:
                self.components = None
                self.data = [zeros((eb.numele,1)) for eb in mesh.eleblx]
            elif self.type == VECTOR:
                self.components = ('x','y','z')[:numdim]
                self.data = [zeros((eb.numele, numdim)) for eb in mesh.eleblx]
            elif self.type == SYMTENSOR:
                self.data = []
                for eb in mesh.eleblx:
                    if numcomp is None:
                        ndir, nshr = eb.eletyp.ndir, eb.eletyp.nshr
                    else:
                        ndir, nshr = {3:(2,1), 4:(3,1), 6:(3,3)}[numcomp]
                    numcomp = ndir + nshr
                    self.components = ('xx','yy','zz')[:ndir]+('xy','yz','xz')[:nshr]
                    self.data.append(zeros((eb.numele, numcomp)))

        else:
            if self.type == SCALAR:
                self.components = None
                self.data = [zeros(eb.numele, eb.npt) for eb in mesh.eleblx]
            elif self.type == VECTOR:
                self.components = ('x','y','z')[:numdim]
                self.data = [zeros((eb.numele, eb.npt, numdim))
                             for eb in mesh.eleblx]
            elif self.type == SYMTENSOR:
                self.data = []
                for eb in mesh.eleblx:
                    if numcomp is None:
                        ndir, nshr = eb.eletyp.ndir, eb.eletyp.nshr
                    else:
                        ndir, nshr = {3:(2,1), 4:(3,1), 6:(3,3)}[numcomp]
                    numcomp = ndir + nshr
                    npt = len(eb.eletyp.gaussp)
                    self.components = ('xx','yy','zz')[:ndir]+('xy','yz','xz')[:nshr]
                    self.data.append(zeros((eb.numele, npt, numcomp)))

    def __getitem__(self, args):
        return self.data[args]

    def add_data(self, data, ix=None, block=None, column=None):
        data = asarray(data)
        if self.position == NODE:
            if ix is None:
                ix = arange(self.data.shape[0])
            if column is not None:
                self.data[ix, column] = data
            else:
                assert data.size == self.data[ix].size, 'Incorrect data shape'
                data.resize(self.data[ix].shape)
                self.data[ix] = data
        else:
            # group element data by block
            if data.ndim == 1:
                data = data.reshape(-1, 1)
            for (i, blk) in enumerate(self.mesh.eleblx):
                ix = [self.mesh.elemap[xel] for xel in blk.labels]
                if column is not None:
                    self.data[i][ix, column] = data[:,0]
                else:
                    self.data[i][ix] = data

    def get_data(self, position=None, block=None):

        if position is None and block is None:
            return array(self.data)

        position = position or self.position

        if block is not None:
            for eb in self.mesh.eleblx:
                if eb.name == block:
                    theblk = eb
                    break
            else:
                raise ValueError('Unrecognized element block {0!r}'.format(block))

        if self.position == ELEMENT:
            if position not in (ELEMENT, ELEMENT_CENTROID):
                raise ValueError('Bad position request for ELEMENT data')
            if block is None:
                return array(self.data)
            return array(self.data[theblk.elems])

        elif self.position == NODE:
            if position != NODE:
                raise ValueError('Bad position request for NODE data')
            if block is not None:
                raise ValueError('Cannot get block data for NODE position')
                ix = sorted(unique(theblk.elecon))
                return self.data[ix]
            return array(self.data)

        # Interpolate Gauss point data to element center
        ave = []
        for (i, dx) in enumerate(self.data):
            et = self.mesh.eleblx[i].eletyp
            ave.append(et.interpolate_to_centroid(dx))
        return ave

    @property
    def values(self):
        if self._values is not None:
            return self._values
        self._values = []
        if self.position == NODE:
            xd = self.data
        else:
            xd = vstack(self.data)
        for (i, row) in enumerate(xd):
            if self.position == NODE:
                label = self.mesh.inodmap[i]
            else:
                label = self.mesh.ielemap[i]
            fv = FieldValue(self.position, label, self.type,
                            self.components, row)
            self._values.append(fv)
        return self._values

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
