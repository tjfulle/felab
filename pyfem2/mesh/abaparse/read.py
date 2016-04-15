import logging
from numpy import *
from collections import OrderedDict
from os.path import basename, isfile

from .parse import AbaqusParser

def mag(a):
    return sqrt(dot(a, a))

def mag2(a):
    return array([mag(x) for x in a])

class Container(OrderedDict):
    def __iter__(self):
        return iter(self.values())
    def get(self, key, default=None):
        K = key.upper()
        for (k, v) in self.items():
            if k.upper() == K:
                return v
        return default

class Node:
    def __init__(self, **kwds):
        self.__dict__.update(**kwds)

class Nodes(Container):
    def __init__(self):
        self._map = {}
        self.labels = []
    def put(self, item):
        nodes = []
        i = len(self)
        for row in item.data:
            label = int(row[0])
            node = Node(label=int(row[0]),
                        coord=array([float(x) for x in row[1:]]))
            self[node.label] = node
            nodes.append(node.label)
            self._map[node.label] = i
            self.labels.append(node.label)
            i += 1
        return nodes

class Element:
    def __init__(self, **kwds):
        self.__dict__.update(**kwds)
    @property
    def centroid(self):
        return average(self.coord, axis=0)
    @property
    def coord(self):
        return array([node.coord for node in self.nodes])
    def normal(self, face):
        et = self.type.upper()
        if et in ('C3D20',):
            assert len(face) == 8
        elif et in ('CAX4',):
            assert len(face) == 2
            xc = array([node.coord for node in self.nodes if node.label in face])
            dx, dy = xc[1,[0,1]] - xc[0,[0,1]]
            n = array([dy, -dx], dtype=float)
            return n / mag(n)
        else:
            raise TypeError('ELEMENT TYPE {0} DOES NOT DEFINE NORMAL'.format(et))

class Elements(Container):
    def __init__(self):
        self._map = {}
        self.labels = []
    def put(self, item, nodes):
        labels = []
        i = len(self)
        et = item.params.get('type').value
        for row in item.data:
            label = int(row[0])
            e_nodes = array([nodes[int(x)] for x in row[1:]])
            e = Element(label=int(row[0]), type=et, nodes=e_nodes,
                        connect=[node.label for node in e_nodes])
            self[e.label] = e
            self._map[e.label] = i
            self.labels.append(e.label)
            labels.append(e.label)
            i += 1
        return labels

class SolidSection:
    def __init__(self, **kwds):
        self.__dict__.update(**kwds)

class SolidSections(Container):
    def put(self, item):
        ss = SolidSection(data=item.data, **item.params.todict())
        self[ss.elset.upper()] = ss
        return ss

class Surface:
    def __init__(self, name):
        self.name = name
        self.faces = []
        self.labels = []
    def put(self, surf):
        for item in surf:
            self.faces.append(item[1])
            self.labels.append(item[0])
    def tolist(self):
        d = dict(zip(('S1','S2','S3','S4','S5','S6','S7','S8'),range(8)))
        surface = []
        for (i, face) in enumerate(self.faces):
            if face.upper().startswith('S'):
                face = eval(face, d, d)
            else:
                face = int(face) - 1
            surface.append((self.labels[i], face))
        return surface

class Surfaces(Container):
    def put(self, item, elsets):
        surf = []
        for row in item.data:
            x, face = row
            labels = elsets.get(x.upper())
            if labels is None:
                labels = [int(x)]
            for label in labels:
                surf.append((label, face))
        p = item.params.get('name')
        name = p.value.strip().upper()
        surface = self.get(name, Surface(name))
        surface.put(surf)
        self[name] = surface
    def todict(self):
        surfaces = {}
        for surface in self.values():
            surfaces[surface.name] = surface.tolist()
        return surfaces

class Material:
    def __init__(self, **kwds):
        self.__dict__.update(**kwds)
        self.defn = {}
    def update(self, item):
        self.defn[item.key.upper()] = item

class Materials(OrderedDict):
    def put(self, item):
        m = Material(**item.params.todict())
        self[m.name.upper()] = m

class NodeSets(Container):
    pass

class ElementSets(Container):
    pass

class AbaqusModel(object):

    def __init__(self, filename):
        self.p = AbaqusParser(lex_optimize=True, yacc_optimize=True)
        self.nodes = Nodes()
        self.elements = Elements()
        self.nodesets = NodeSets()
        self.elsets = ElementSets()
        self.surfaces = Surfaces()
        self.solid_sections = SolidSections()
        self.materials = Materials()
        self._c = None

        self.parse(filename)

    @property
    def elemap(self):
        return self.elements._map

    @property
    def nodmap(self):
        return self.nodes._map

    def get_nodes_in_set(self, nset, disp=0):
        labels = self.nodesets.get(nset.upper())
        if disp:
            return [self.nodes[label] for label in labels]
        return labels

    def get_elements_in_set(self, elset, disp=0):
        labels = self.elsets.get(elset.upper())
        if disp:
            return [self.elements[label] for label in labels]
        return labels

    def elem_centroids(self, labels=None):
        if self._c is None:
            self._c = array([self.elements[label].centroid
                             for label in self.elements.labels])
        if labels is None:
            return self._c
        ix = [self.elemap[label] for label in labels]
        return self._c[ix]

    def element_table(self):
        return [[el.label]+el.connect for el in self.elements]

    def node_table(self):
        return [[node.label]+node.coord.tolist() for node in self.nodes]

    def element_blocks(self, fun=None):
        """Form ExodusII type element blocks"""
        # CHECK SOLID SECTIONS FOR CONSISTENCY WITH ELEMENT BLOCK REQUIREMENT
        # OF SINGLE ELEMENT TYPE
        eleblx = {}
        for ss in self.solid_sections:
            labels = self.elsets.get(ss.elset)
            elems = [self.elements[label] for label in labels]
            eletyp = set([e.type for e in elems])
            if len(eletyp) != 1:
                raise ValueError('ELEMENT BLOCKS FORMED FROM SOLID SECTIONS '
                                 'MUST CONTAIN ONLY ONE ELEMENT TYPE')
            et = elems[0].type if fun is None else fun(elems[0].type)
            eleblx[ss.elset] = (et, labels)
        return eleblx

    def get_closest_element(self, pos, disp=0):
        xc = self.elem_centroids()
        i = argmin(mag2(abs(xc - pos)))
        label = self.elements.labels[i]
        if disp:
            return self.elements[label]
        return label

    def shift(self, offset):
        for (label, node) in self.nodes.items():
            node.coord += offset
        self._c = None

    def parse(self, filename):
        assert isfile(filename)
        self.f = basename(filename)
        buf = open(filename).read()
        t = self.p.parse(buf, self.f, debuglevel=0)

        notread = []
        for item in t:

            if item.key == 'node':
                labels = self.nodes.put(item)
                p = item.params.get('nset')
                if p is not None:
                    self.nodesets.setdefault(p.value.upper(), []).extend(labels)

            elif item.key == 'element':
                labels = self.elements.put(item, self.nodes)
                p = item.params.get('elset')
                if p is not None:
                    self.elsets.setdefault(p.value.upper(), []).extend(labels)

            elif item.key == 'solid_section':
                ss = self.solid_sections.put(item)
                for label in self.elsets[ss.elset.upper()]:
                    self.elements[label].solid_section = ss

            elif item.key == 'elset':
                name, labels = self.parse_set(item)
                self.elsets.setdefault(name.upper(), []).extend(labels)

            elif item.key == 'nset':
                name, labels = self.parse_set(item)
                self.nodesets.setdefault(name.upper(), []).extend(labels)

            elif item.key == 'surface':
                self.surfaces.put(item, self.elsets)

            elif item.key == 'material':
                self.materials.put(item)

            elif item.key in ('elastic', 'density', 'expansion'):
                mat = list(self.materials.values())[-1]
                mat.update(item)

            else:
                notread.append(item.key)

        if notread:
            logging.debug('THE FOLLOWING KEYWORDS AND THEIR DATA WERE NOT READ:\n'
                          '{0}'.format(', '.join(notread)))

    def parse_set(self, item):
        generate = 'generate' in item.params
        if not generate:
            labels = [int(e) for row in item.data for e in row]
        else:
            labels = []
            for row in item.data:
                start, stop, step = [int(n) for n in row]
                labels.extend(range(start, stop+1, step))
        p = item.params.get(item.key)
        return p.value, labels


if __name__ == '__main__':
    filename = 'mmxmn.inp'
    assert isfile(filename)
    a = AbaqusInput(filename)
    a.shift([-1., 0., 0.])
    el = a.elements[1]
    xc = el.centroid
    el1 = a.get_closest_element(xc, disp=1)
    print(el1.label)
