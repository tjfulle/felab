import re
import os
import sys
import logging
import datetime
from numpy import *

SCALAR = 'Scalar'
VECTOR = 'Vector'
SYMTENSOR = 'Symmetric Tensor'
INTEGRATION_POINT = 'Integration Point'
ELEMENT_CENTROID = 'Element Centroid'

class FiniteElementDatabase(object):

    def __init__(self, num_dim, node_map, node_variables, elem_map, elem_blx):

        self.num_elem = len(elem_map)
        self.num_node = len(node_map)

        # Node map: nodmap1[i] is the external node label of the ith node
        self.node_map = node_map
        self.elem_map = elem_map

        self.count = 0
        self.time = zeros((2, 3))
        self.time_step = zeros(2)

        self._initialize_node_data(node_variables)
        self._initialize_elem_data(elem_blx)

    def _initialize_node_data(self, node_variables):
        # ------------------------------------------------------------------- #
        # ----------------------------------------------------- node data --- #
        # ------------------------------------------------------------------- #
        m, data = 0, []
        self.node_var_columns = {}
        self.node_var_names = []
        for node_variable in node_variables:
            name, type = node_variable[:2]
            a = zeros(self.num_node)
            if len(node_variable) > 2:
                # Initial value
                a[:] = node_variable[2]
            if type == SCALAR:
                data.append(a)
                self.node_var_names.append(name)
                self.node_var_columns[name] = (m, m+1)
                m += 1
            elif type == VECTOR:
                data.extend([a] * num_dim)
                keys = self.vector_components(name, num_dim)
                self.node_var_names.extend(keys)
                self.node_var_columns[name] = (m, len(keys))
                for key in keys:
                    self.node_var_columns[key] = (m, m+1)
                    m += 1
            else:
                raise TypeError('Unknown node variable type')
        data = column_stack(data)
        self.node_data = repeat(data, 2, axis=0).reshape(2, self.num_node, -1)

    def _initialize_elem_data(self, elem_blx):

        # ------------------------------------------------------------------- #
        # -------------------------------------------------- element data --- #
        # ------------------------------------------------------------------- #
        self.elem_blk_data = {}
        self.elem_var_columns = {}
        self.elem_var_names = {}
        for eb in elem_blx:
            m, data = 0, []
            self.elem_blk_data[eb.name] = {}
            self.elem_var_columns[eb.name] = {}
            self.elem_var_names[eb.name] = []

            if not eb.eletyp.variables():
                continue

            num_elem_in_blk = len(eb.labels)
            ndir, nshr = eb.eletyp.ndir, eb.eletyp.nshr
            ntens = ndir + nshr

            if eb.eletyp.integration:
                position = INTEGRATION_POINT
            else:
                position = ELEMENT_CENTROID

            for elem_variable in eb.eletyp.variables():
                name, type = elem_variable[:2]
                idata = 0. if len(elem_variable) <= 2 else elem_variable[2]

                if type == SCALAR:
                    a = zeros(num_elem_in_blk)
                    a[:] = idata
                    data.append(a)
                    self.elem_var_names[eb.name].append(name)
                    self.elem_var_columns[eb.name][name] = (m, m+1)
                    m += 1
                elif type == VECTOR:
                    a = zeros(num_node)
                    a[:] = idata
                    data.extend([a] * num_dim)
                    keys = self.vector_components(name, num_dim)
                    self.elem_var_names[eb.name].extend(keys)
                    self.elem_var_columns[eb.name][name] = (m, len(keys))
                    for key in keys:
                        self.elem_var_columns[eb.name][key] = (m, m+1)
                        m += 1
                elif type == SYMTENSOR:
                    a = zeros(ntens)
                    if abs(idata - 1.) < 1e-12:
                        idata = [1.] * ndir + [0.] * nshr
                    a[:] = idata
                    a = row_stack([a] * num_elem_in_blk)
                    keys = self.symtensor_components(name, ndir, nshr)
                    self.elem_var_names[eb.name].extend(keys)
                    self.elem_var_columns[name] = (m, len(keys))
                    for key in keys:
                        self.elem_var_columns[eb.name][key] = (m, m+1)
                        m += 1
                    data.extend([x for x in a.T])
                else:
                    raise TypeError('Unknown element variable type')

            data = column_stack(data)
            if not eb.eletyp.integration:
                shape = (2, num_elem_in_blk, -1)
                eb_data = repeat(data, 2, axis=0).reshape(shape)
            else:
                ne, np = num_elem_in_blk, eb.eletyp.integration
                eb_data = repeat(data, np, axis=0).reshape(ne, np, -1)
                eb_data = repeat(eb_data, 2, axis=0).reshape(2,ne,np,-1)

            self.elem_blk_data[eb.name] = eb_data

        return

    def get_node_data(self, name=None, node=None, label=None):
        if name is not None:
            a, b = self.node_var_columns[name]
            return self.node_data[0, :, a:b]
        if node is not None:
            i = node
        elif label is not None:
            i = self.node_map[label]
        return self.node_data[0, i, :]

    def get_elem_data(self, block, ip, name=None, element=None, label=None):
        d = self.elem_blk_data[block]
        if name is not None:
            a, b = self.elem_var_columns[name]
            return d[0, :, ip, a:b]
        if element is not None:
            i = element
        elif label is not None:
            i = self.elem_map[label]
        return d[0, i, ip, :]

    def advance(self):
        self.node_data[0] = self.node_data[1]
        for (eb, eb_data) in self.elem_blk_data.items():
            eb_data[0] = eb_data[1]

    def vector_components(self, name, n):
        return ['{0}_{1}'.format(name, x) for x in 'XYZ'[:n]]

    def symtensor_components(self, name, ndir, nshr):
        keys  = ['{0}_{1}'.format(name, x) for x in ['XX', 'YY', 'ZZ'][:ndir]]
        keys += ['{0}_{1}'.format(name, x) for x in ['XY', 'YZ', 'XZ'][:nshr]]
        return keys

if __name__ == '__main__':
    class Temp:
        pass

    num_dim = 3
    node_map = {1: 0, 2: 1, 3: 2}
    elem_map = {1: 0, 2: 1, 3: 2}
    elem_blx = []
    node_variables = (('U', VECTOR), ('T', SCALAR, 295))

    def variables():
        return (('S', SYMTENSOR), ('V', SYMTENSOR, 1), ('H', SCALAR, 2))
    eb1 = Temp()
    eb1.name = 'B1'
    eb1.labels = [1, 2, 3]
    eb1.eletyp = Temp()
    eb1.eletyp.ndir = 2
    eb1.eletyp.nshr = 1
    eb1.eletyp.integration = 2
    eb1.eletyp.variables = variables

    eb2 = Temp()
    eb2.name = 'B2'
    eb2.labels = [1, 2, 3]
    eb2.eletyp = Temp()
    eb2.eletyp.ndir = 2
    eb2.eletyp.nshr = 1
    eb2.eletyp.integration = 2
    eb2.eletyp.variables = variables

    elem_blx = [eb1, eb2]

    dw = FiniteElementDatabase(num_dim, node_map, node_variables,
                               elem_map, elem_blx)


    print(dw.node_data)
    u = dw.get_node_data('U')
    print(u)
    u[:,0] += 1
    u[:,1] += 2
    u[:,2] += 3
    u1 = dw.get_node_data('U')
    print(u1)
    print(dw.node_var_columns)
    u2 = dw.get_node_data('U_Y')
    print(u2)

    print(dw.get_node_data('T'))
    print(dw.get_node_data(node=0))
    print(dw.get_node_data(label=3))

    dw.advance()
    for a,b in dw.elem_blk_data.items():
        print(a,b)
        print(b.shape)
