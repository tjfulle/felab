import logging
import pickle
from numpy import *
from numpy.ma import masked, masked_array
from collections import OrderedDict
from os.path import basename, isfile, splitext

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
        """Find the normal to the plane defined by points"""
        et = self.type.upper()
        xc = array([node.coord for node in self.nodes if node.label in face])
        if et[:5] == 'C3D20':
            assert len(xc) == 8
            x1 = xc[1] - xc[0]
            for x in xc[2:]:
                x2 = x - xc[0]
                if abs(dot(x1, x2)) < 1e-4:
                    continue
                n = cross(x2, x1)
                return n / mag(n)
        elif et[:4] == 'C3D8':
            assert len(face) == 4
            x1 = xc[1] - xc[0]
            for x in xc[1:]:
                x2 = x - xc[0]
                if abs(dot(x1, x2)) < 1e-4:
                    continue
                n = cross(x2, x1)
                return n / mag(n)
        elif et in ('CAX4', 'PE4', 'PS4'):
            assert len(face) == 2
            dx, dy = xc[1,[0,1]] - xc[0,[0,1]]
            n = array([dy, -dx], dtype=float)
            return n / mag(n)
        else:
            raise TypeError('ELEMENT TYPE {0} DOES NOT DEFINE NORMAL'.format(et))
        raise ValueError('UNABLE TO FIND NORMAL')
    def radius(self, axis):
        ij = [x for x in range(3) if x != axis]
        return mag(self.centroid[ij])

class Elements(Container):
    def __init__(self):
        self._map = {}
        self.labels = []
    def put(self, item, nodes):
        labels = []
        i = len(self)
        et = item.params.get('type')
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
        name = item.params.get('name').upper()
        surface = self.get(name, Surface(name))
        surface.put(surf)
        self[name] = surface
    def todict(self):
        surfaces = {}
        for surface in self.values():
            surfaces[surface.name] = surface.tolist()
        return surfaces

class Material:
    def __init__(self, item):
        self.card = item
        self.name = item.params.get('name').upper()
        self.subcards = []
    def update(self, item):
        self.subcards.append(item)
    def toinp(self):
        o = self.card.toinp()
        for card in self.subcards:
            o += '\n' + card.toinp()
        return o.strip()

class Materials(OrderedDict):
    def put(self, item):
        m = Material(item)
        self[m.name] = m

class Orientations(OrderedDict):
    def put(self, item):
        pass

class SetContainer(Container):
    def put(self, item, labels=None):
        if labels is not None:
            self.setdefault(item.upper(), []).extend(labels)
            return

        generate = 'generate' in item.params
        if not generate:
            try:
                labels = [int(e) for row in item.data for e in row]
            except ValueError:
                xsets1 = [str(e).upper() for row in item.data for e in row]
                labels = [label for xset in xsets1 for label in self[xset]]
        else:
            labels = []
            for row in item.data:
                start, stop, step = [int(n) for n in row]
                labels.extend(range(start, stop+1, step))
        labels = sorted(list(set(labels)))
        name = item.params.get(item.key).upper()
        self.setdefault(name, []).extend(labels)

class NodeSets(SetContainer):
    pass

class ElementSets(SetContainer):
    pass

def AbaqusModelFactory(filename):
    root, ext = splitext(filename)
    pkl = root + '.p'
    if isfile(pkl):
        with open(pkl, 'rb') as fh:
            model = pickle.load(fh)
    else:
        model = AbaqusModel(filename)
        with open(pkl, 'wb') as fh:
            pickle.dump(model, fh)
    return model

class AbaqusModel(object):

    def __init__(self, filename):

        self.nodes = Nodes()
        self.elements = Elements()
        self.nodesets = NodeSets()
        self.elsets = ElementSets()
        self.surfaces = Surfaces()
        self.orientations = Orientations()
        self.solid_sections = SolidSections()
        self.materials = Materials()
        self._c = None
        self._x = None

        self.parse(filename)

    @property
    def elemap(self):
        return self.elements._map

    @property
    def nodmap(self):
        return self.nodes._map

    @property
    def dimension(self):
        return self.get_coord().shape[1]

    def get_min_on_axis(self, axis):
        return amin(self.get_coord()[:,axis])

    def get_max_on_axis(self, axis):
        return amax(self.get_coord()[:,axis])

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

    def get_elem_centroids(self, labels=None, subset=None):
        if self._c is None:
            self._c = array([self.elements[label].centroid
                             for label in self.elements.labels])
        if labels is not None:
            ix = [self.elemap[label] for label in labels]
            return self._c[ix]
        if subset is not None:
            ix = [self.elemap[e.label] for e in subset]
            return self._c[ix]
        return self._c

    def get_coord(self, labels=None):
        if self._x is None:
            self._x = array([self.nodes[label].coord
                             for label in self.nodes.labels])
        if labels is None:
            return self._x
        ix = [self.nodmap[label] for label in labels]
        return self._x[ix]

    def get_length(self, axis):
        xc = self.get_elem_centroids()[:,axis]
        return abs(amax(xc) - amin(xc))

    def get_radius(self, axis, pos=None, subset=None, at_min=0, at_max=0):
        xc = self.get_elem_centroids()
        if subset is not None:
            ix = [self.elemap[e.label] for e in subset]
            xc = xc[ix]
            if at_min:
                i = argmin(xc[:,axis])
            elif at_max:
                i = argmax(xc[:,axis])
            else:
                raise TypeError('Requires at_min or at_max keyword')
        else:
            if pos is None:
                raise TypeError('pos keyword required')
            i = argmin(abs(xc[:,axis]-pos))

        r = sqrt(sum([xc[i, j] ** 2 for j in range(xc.shape[1]) if j!=axis]))
        return r

    def get_closest_element(self, pos, disp=0, axis=0, radial=0,
                            skip_sets=None, skip_labels=None):

        skip_labels = skip_labels or []
        skip_sets = skip_sets or []
        for xset in enumerate(skip_sets):
            skip_labels.extend(self.elsets[xset.upper()])

        xc = self.get_elem_centroids()
        clen = abs(average(diff(sorted(xc[:,axis]))))

        if radial:
            # RADIAL COORDINATES
            if len(pos) == 1:
                raise ValueError('POS MUST BE AT LEAST 2 DIMENSIONS')

            z = pos[axis]
            if len(pos) == 2:
                r = pos[{0:1, 1:0}[axis]]
            else:
                r = sqrt(sum(pos[i]**2 for i in range(3) if i!=axis))
            pos = array([r, z])

            # ALL ELEMENTS CENTROIDS
            z = xc[:,axis]
            if xc.shape[1] == 2:
                r = xc[:, {0:1, 1:0}[axis]]
            else:
                r = sqrt(sum(xc[:,i]**2 for i in range(3) if i!=axis))
            xc = column_stack((r, z))

        dx = mag2(abs(xc - pos))

        if skip_labels:
            skip_labels = sorted(list(set(skip_labels)))
            ix = [self.elemap[label] for label in skip_labels]
            dx_m = masked_array(dx)
            dx_m[ix] = masked
            index = dx_m.argmin()

        else:
            index = dx.argmin()

        label = self.elements.labels[index]
        if disp:
            return self.elements[label]
        return label

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

    def shift(self, offset):
        for (label, node) in self.nodes.items():
            node.coord += offset
        self._c = None
        self._x = None

    def parse(self, filename):
        assert isfile(filename)
        self.f = basename(filename)
        p = AbaqusParser(lex_optimize=False, yacc_optimize=False)
        buf = open(filename).read()
        t = p.parse(buf, self.f, debuglevel=0)

        notread = []
        for item in t:

            if item.key == 'node':
                labels = self.nodes.put(item)
                p = item.params.get('nset')
                if p is not None:
                    self.nodesets.put(p.upper(), labels=labels)

            elif item.key == 'element':
                labels = self.elements.put(item, self.nodes)
                p = item.params.get('elset')
                if p is not None:
                    self.elsets.put(p.upper(), labels=labels)

            elif item.key == 'solid_section':
                ss = self.solid_sections.put(item)
                for label in self.elsets[ss.elset.upper()]:
                    self.elements[label].solid_section = ss

            elif item.key == 'elset':
                self.elsets.put(item)

            elif item.key == 'nset':
                elset = item.params.get('elset')
                self.nodesets.put(item)

            elif item.key == 'surface':
                self.surfaces.put(item, self.elsets)

            elif item.key == 'material':
                self.materials.put(item)

            elif item.key in ('hyperelastic', 'elastic', 'density', 'expansion'):
                mat = list(self.materials.values())[-1]
                mat.update(item)

            elif item.key == 'orientation':
                self.orientations.put(item)

            elif item.key == 'heading':
                self.heading = item.data

            else:
                notread.append(item.key)

        if notread:
            logging.debug('THE FOLLOWING KEYWORDS AND THEIR DATA WERE NOT READ:\n'
                          '{0}'.format(', '.join(notread)))

def test():
    from os.path import dirname, realpath, join
    D = dirname(realpath(__file__))
    filename = join(D, 'mmxmn.inp')
    assert isfile(filename)
    a = AbaqusModel(filename)
    a.shift([-1., 0., 0.])
    el = a.elements[1]
    xc = el.centroid
    el1 = a.get_closest_element(xc, disp=1)
    print(el1.label)
    print(list(a.materials.values())[0].toinp())
