#!/usr/bin/env python
import datetime
import os
from os.path import join, isdir
from numpy import *
import xml.dom.minidom as xdom

from .constants import *
from .elemlib1 import ElementFamily

# Linear cells
VTK_EMPTY_CELL       = 0
VTK_VERTEX           = 1
VTK_POLY_VERTEX      = 2
VTK_LINE             = 3
VTK_POLY_LINE        = 4
VTK_TRIANGLE         = 5
VTK_TRIANGLE_STRIP   = 6
VTK_POLYGON          = 7
VTK_PIXEL            = 8
VTK_QUAD             = 9
VTK_TETRA            = 10
VTK_VOXEL            = 11
VTK_HEXAHEDRON       = 12
VTK_WEDGE            = 13
VTK_PYRAMID          = 14
VTK_PENTAGONAL_PRISM = 15
VTK_HEXAGONAL_PRISM  = 16

# Quadratic, isoparametric cells
VTK_QUADRATIC_EDGE                   = 21
VTK_QUADRATIC_TRIANGLE               = 22
VTK_QUADRATIC_QUAD                   = 23
VTK_QUADRATIC_TETRA                  = 24
VTK_QUADRATIC_HEXAHEDRON             = 25
VTK_QUADRATIC_WEDGE                  = 26
VTK_QUADRATIC_PYRAMID                = 27
VTK_BIQUADRATIC_QUAD                 = 28
VTK_TRIQUADRATIC_HEXAHEDRON          = 29
VTK_QUADRATIC_LINEAR_QUAD            = 30
VTK_QUADRATIC_LINEAR_WEDGE           = 31
VTK_BIQUADRATIC_QUADRATIC_WEDGE      = 32
VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33

VTK_S = 'VTK'
VTK_USGRID = 'UnstructuredGrid'
VTK_DATARR = 'DataArray'
VTK_FILE = 'VTKFile'

INDENT = 0
def arrtostr2(a, fmt='.18f', indent='', newl='\n'):
    a1 = ['{0}{1}'.format(indent, arrtostr(row, fmt=fmt, newl='')) for row in a]
    return '{0}{1}{0}{2}'.format(newl, newl.join(a1), indent[:-INDENT])

def arrtostr(a, fmt='.18f', newl='\n', indent=''):
    a1 = ' '.join('{0:{1}}'.format(x, fmt) for x in a)
    return '{0}{3}{1}{0}{2}'.format(newl, a1, indent[:-INDENT], indent)

def vtk_cell_types(et=None, dimensions=None, nodes=None):
    if et is not None:
        dimensions, nodes = et.dimensions, et.nodes
    if nodes == 2:
        return VTK_LINE
    if dimensions == 1:
        return VTK_LINE
    elif dimensions == 2:
        if nodes == 3:
            return VTK_TRIANGLE
        elif nodes == 4:
            return VTK_QUAD
        elif nodes == 6:
            return VTK_QUADRATIC_TRIANGLE
        elif nodes == 8:
            return VTK_QUADRATIC_QUAD

def pyfem_elem_type(et):
    dimensions, nodes = {VTK_LINE: (None, 2),
                         VTK_TRIANGLE: (2, 3),
                         VTK_QUAD: (2, 4),
                         VTK_QUADRATIC_QUAD: (2, 8)}[et]
    return ElementFamily(dimensions, nodes)

def arrtotens(a):
    if a.ndim == 1:
        if a.size == 4:
            return array((a[0],a[3],0,a[3],a[1],0,0,0,a[2]))
        elif a.size == 6:
            return array((a[0],a[3],a[5],a[3],a[1],a[4],a[5],a[4],a[2]))
    if a.shape[1] == 4:
        z = zeros(a.shape[0])
        return column_stack((a[:,0],a[:,3],z,
                             a[:,3],a[:,1],z,
                             z,     z,     a[:,2]))
    elif a.shape[1] == 6:
        return column_stack((a[:,0],a[:,3],a[:,5],
                             a[:,3],a[:,1],a[:,4],
                             a[:,5],a[:,4],a[:,2]))
    raise ValueError('Unknown array size')

class VTKFile(object):
    def __init__(self, jobid=None, filename=None, mode='r'):
        self.jobid = jobid
        self.mode = mode.lower()
        self.pvd = None
        self.count = 0
        if mode not in 'rw':
            raise ValueError('Unknown file mode')
        if self.mode[0] == 'w':
            self.datadir = self.jobid + '.vtu.d'
            self.pvd = xdom.Document()
            self.pvd_root = self.pvd.createElementNS('VTK', 'VTKFile')
            self.pvd_root.setAttribute('type', 'Collection')
            self.pvd_root.setAttribute('version', '0.1')
            self.pvd_root.setAttribute('byte_order', 'LittleEndian')
            self.pvd.appendChild(self.pvd_root)
            self.pvd_coll = self.pvd.createElementNS('VTK', 'Collection')
            self.pvd_root.appendChild(self.pvd_coll)
        else:
            if filename is None:
                raise ValueError('No filename given')
            if not os.path.isfile(filename):
                f = os.path.join(DATA_D, filename)
                if not os.path.isfile(f):
                    raise ValueError('File not found: {0!r}'.format(filename))
                filename = f
            self.filename = filename
            x = self.readmesh(filename)
            (self.coord, self.nodlab, self.elecon, self.elelab, self.nodesets,
             self.elemsets, self.surfaces) = x

    def put_init(self, coord, nodlab, elelab, eletyp, elecon):
        self.coord = asarray(coord)
        self.nodlab = asarray(nodlab, dtype=int)
        self.elecon = asarray(elecon, dtype=int)
        self.elelab = asarray(elelab, dtype=int)
        self.eletyp = asarray(eletyp, dtype=object)
        self.nodes, self.dimensions = coord.shape
        self.numele, self.max_num_node_per_elem = elecon.shape

    def readmesh(self, filename):
        if filename.endswith('.gz'):
            f = gzip.open(filename, 'rb')
            self.doc = xdom.parseString(f.read())
            f.close()
        else:
            self.doc = xdom.parse(filename)
        root = self.doc.getElementsByTagName('VTKFile')[0]
        grid = root.getElementsByTagName('UnstructuredGrid')
        if not grid:
            logging.warn('NO GRID INFORMATION FOUND in {0!r}'.format(filename))
            return None
        if len(grid) > 1:
            logging.warn('MULTIPLE GRIDS NOT SUPPORTED')
            return None
        grid = grid[0]
        # Piece 0 (only one)
        piece = grid.getElementsByTagName('Piece')[0]
        nodes = int(piece.getAttribute('NumberOfPoints'))
        numele = int(piece.getAttribute('NumberOfCells'))
        points = piece.getElementsByTagName('Points')[0]
        el = points.getElementsByTagName('DataArray')[0]
        numcmp = int(el.getAttribute('NumberOfComponents'))
        s = [line for line in el.firstChild.data.split('\n') if line.split()]
        coord = array([float(a) for o in s for a in o.split()]).reshape(nodes, -1)
        assert coord.shape[1] == numcmp, 'Error reading coord shape'

        try:
            dimensions = int(el.getAttribute('NumberOfProblemDimensions'))
        except ValueError:
            dimensions = 0
        if not dimensions:
            dimensions = coord.shape[1]
            # check dimension
            if dimensions == 3 and all(abs(coord[:,dimensions-1] < 1e-12)):
                dimensions = 2
        coord = coord[:, :dimensions]

        # Cells
        cells = piece.getElementsByTagName('Cells')[0]
        data = cells.getElementsByTagName('DataArray')
        for item in data:
            s = lambda el: [a for x in el.firstChild.data.split('\n')
                            for a in x.split() if x.split()]
            name = item.getAttribute('Name')
            if name == 'connectivity':
                conn = array([int(x) for x in s(item)])
            elif name == 'offsets':
                offsets = array([int(x) for x in s(item)])
            elif name == 'types':
                eletyp = array([pyfem_elem_type(int(x)) for x in s(item)])
        if len(offsets) > 1:
            maxnod = max(diff(offsets))
        else:
            maxnod = offsets[0]

        # format connectivity
        elecon = zeros((numele,maxnod), dtype=int)
        k = 0
        for (el, offset) in enumerate(offsets):
            j = offset - k
            elecon[el,:j] = [n for n in conn[k:k+j]]
            k += j

        # Field data
        fd = root.getElementsByTagName('FieldData')

        elemsets, surfaces, nodesets = {}, {}, {}
        nodlab, elelab = range(nodes), range(numele)

        if not fd:
            return (coord, nodlab, elecon, elelab,
                    elemsets, surfaces, nodesets)
        fd = fd[0]

        # Node table
        node_labels = fd.getElementsByTagName('NodeLabels')
        if node_labels:
            el = node_labels[0].getElementsByTagName('DataArray')[0]
            nodlab = [int(a) for x in el.firstChild.data.split('\n')
                      for a in x.split() if x.split()]

        # Element table
        elem_labels = fd.getElementsByTagName('ElementLabels')
        if elem_labels:
            el = elem_labels[0].getElementsByTagName('DataArray')[0]
            elelab = [int(a) for x in el.firstChild.data.split('\n')
                      for a in x.split() if x.split()]

        # Element sets
        for item in fd.getElementsByTagName('ElementSet'):
            name = item.getAttribute('Name')
            el = item.getElementsByTagName('DataArray')[0]
            if el.getAttribute('type') == 'String': f = str
            else: f = int
            elems = [f(a) for x in el.firstChild.data.split('\n')
                     for a in x.split() if x.split()]
            elemsets[name] = elems

        # Surfaces
        multisurf = []
        for item in fd.getElementsByTagName('Surface'):
            name = item.getAttribute('Name')
            el = item.getElementsByTagName('DataArray')[0]
            lines = [s.split() for s in el.firstChild.data.split('\n') if s.split()]
            if el.getAttribute('type') == 'String':
                multisurf.append(name)
                surf = [str(a) for line in lines for a in line]
            else:
                surf = [[int(a) for a in line] for line in lines]
            surfaces[name] = [x for x in surf if x]
        for name in multisurf:
            surfs = []
            for surf in surfaces[name]:
                surfs.extend(surfaces[surf])
            surfaces[name] = surfs

        # Node sets
        for item in fd.getElementsByTagName('NodeSet'):
            name = item.getAttribute('Name')
            el = item.getElementsByTagName('DataArray')[0]
            if el.getAttribute('type') == 'String': f = str
            else: f = int
            nodes = [f(a) for x in el.firstChild.data.split('\n')
                     for a in x.split() if x.split()]
            nodesets[name] = nodes

        return coord, nodlab, elecon, elelab, nodesets, elemsets, surfaces

    def snapshot(self, u=None, **kwds):
        if not isdir(self.datadir):
            os.makedirs(self.datadir)
        filename = join(self.datadir, self.jobid + '-{0:04d}.vtu'.format(self.count))
        self.write_vtu_file(filename=filename, u=u, **kwds)
        # Write the updated pvd file
        ds = self.pvd.createElementNS('VTK', 'DataSet')
        ds.setAttribute('timestep', str(self.count))
        ds.setAttribute('group', '')
        ds.setAttribute('part', '0')
        ds.setAttribute('file', filename)
        self.pvd_coll.appendChild(ds)
        with open(self.jobid + '.pvd', 'w') as fh:
            fh.write(self.pvd.toprettyxml(indent=''))
        self.count += 1

    def write_vtu_file(self, filename=None, u=None, **kwds):
        filename = filename or self.jobid + '.vtu'
        self.doc = xdom.Document()
        root = self.doc.createElementNS('VTK', 'VTKFile')
        root.setAttribute('type', 'UnstructuredGrid')
        root.setAttribute('version', '0.1')
        root.setAttribute('byte_order', 'LittleEndian')
        self.doc.appendChild(root)

        # unstructured grid
        self.write_grid(u=u)
        if u is not None:
            kwds['u'] = u

        # Data at nodes
        pd = self.create_element('PointData', parent='Piece')
        cd = self.create_element('CellData', parent='Piece')

        for (key, val) in kwds.items():
            val = asarray(val)
            da = self.doc.createElementNS('VTK', 'DataArray')
            da.setAttribute('Name', key)
            if val.ndim == 1:
                # Scalar data
                nc = 1
                val = val.reshape(-1, 1)
            else:
                if val.shape[1] == self.dimensions:
                    nc = 3
                    # Vector data
                    if val.shape[1] != 3:
                        z = zeros((self.nodes,3-val.shape[1]))
                        val = column_stack((val, z))
                elif val.shape[1] == 4:
                    val = arrtotens(val, 1)
                    nc = 6
                else:
                    nc = val.shape[1]

            da.setAttribute('NumberOfComponents', str(nc))
            da.setAttribute('type', 'Float32')
            da.setAttribute('format', 'ascii')
            if len(val) == self.nodes:
                pd.appendChild(da)
            elif len(val) == self.numele:
                cd.appendChild(da)
            else:
                raise Exception('unknown data shape')
            da.appendChild(self.doc.createTextNode(arrtostr2(val)))

        with open(filename, 'w') as fh:
            fh.write(self.doc.toprettyxml(indent=''))

    def flatten(self, the_el=None):
        if the_el is None:
            the_el = self.doc
        flattened = []
        for el in the_el.childNodes:
            if el.nodeType != 1:
                continue
            flattened.append(el)
            flattened.extend(self.flatten(el))
        return flattened

    def create_element(self, name, parent='VTKFile', **kwds):
        el = self.doc.createElementNS('VTK', name)
        for (key, val) in kwds.items():
            el.setAttribute(key, str(val))
        for el1 in self.flatten():
            if el1.nodeName == parent:
                el1.appendChild(el)
                break
        else:
            raise ValueError('Parent node not found')
        return el

    def get_element(self, name, parent=None):
        for el in self.flatten():
            if el.nodeName == name:
                return el
        if parent is not None:
            return self.create_element(name, parent=parent)
        raise ValueError('No element {0!r}'.format(name))

    def write_grid(self, u=None):

        # write the unstructured grid

        # Unstructured grid element
        usgrid = self.create_element('UnstructuredGrid')

        # Piece 0 (only one)
        piece = self.create_element('Piece', parent=usgrid.nodeName,
                                    NumberOfPoints=self.nodes,
                                    NumberOfCells=self.numele)

        # Points
        points = self.create_element('Points', parent=piece.nodeName)

        # Point location data
        da = self.create_element('DataArray', parent=points.nodeName,
                                 type='Float32', format='ascii',
                                 NumberOfComponents=3,
                                 NumberOfProblemDimensions=self.dimensions)
        x = array(self.coord)
        if u is not None:
            x += u
        if self.dimensions != 3:
            z = zeros(self.nodes)
            for i in range(3-self.dimensions):
                x = column_stack((x, z))
        da.appendChild(self.doc.createTextNode(arrtostr2(x, indent='')))

        # Cells
        cells = self.create_element('Cells', parent=piece.nodeName)

        # Cell connectivity
        da = self.create_element('DataArray', parent=cells.nodeName,
                                 type='Int32', Name='connectivity', format='ascii')
        o = [[n for n in el[:self.eletyp[e].nodes]]
             for (e, el) in enumerate(self.elecon)]
        da.appendChild(self.doc.createTextNode(arrtostr2(o, 'd', indent='')))

        # Cell offsets
        da = self.create_element('DataArray', parent=cells.nodeName,
                                 type='Int32', Name='offsets', format='ascii')
        o, k = [], 0
        for (e, el) in enumerate(self.elecon):
            k += self.eletyp[e].nodes
            o.append(k)
        da.appendChild(self.doc.createTextNode(arrtostr(o, 'd', indent='')))

        # Cell types
        da = self.create_element('DataArray', parent=cells.nodeName,
                                 type='UInt8', Name='types', format='ascii')
        o = [vtk_cell_types(et) for et in self.eletyp]
        da.appendChild(self.doc.createTextNode(arrtostr(o, 'd', indent='')))

        # Node and element tables
        self.create_fd('NodeLabels', 'Int32', self.nodlab)
        self.create_fd('ElementLabels', 'Int32', self.elelab)

        return

    def create_fd(self, name, dtype, arr):

        # Field data
        fd = self.get_element('FieldData', parent='VTKFile')

        el = self.create_element(name, parent=fd.nodeName)
        da = self.create_element('DataArray', parent=el.nodeName,
                                 type=dtype, format='ascii')
        arr = asarray(arr)
        fmt = {'Int32': 'd', 'Float32': '.18f'}[dtype]
        if arr.ndim == 1:
            string = arrtostr(arr, fmt, indent='')
        else:
            string = arrtostr2(arr, fmt, indent='')
        da.appendChild(self.doc.createTextNode(string))
        return

    def close(self):
        self.pvd.close()

def WriteVTUMesh(filename, coord, nodlab, elelab, eletyp, elecon, check=0):
    jobid = os.path.splitext(filename)[0]
    f = VTKFile(jobid, mode='w')
    f.put_init(coord, nodlab, elelab, eletyp, elecon)
    f.write_vtu_file()
    if check:
        x, nl, ec, el = ReadMesh(filename, disp=1)
        assert allclose(x, coord), 'coord'
        assert allclose(nl, nodlab), 'nodlab'
        assert allclose(ec, elecon), 'elecon'
        assert allclose(el, elelab), 'elelab'

def WriteFEResults(jobid, coord, nodmap, elemap, eletyp, elecon, u=None, **kwds):
    nodlab = sorted(nodmap.keys(), key=lambda k: nodmap[k])
    elelab = sorted(elemap.keys(), key=lambda k: elemap[k])
    f = VTKFile(jobid, mode='w')
    f.put_init(coord, nodlab, elelab, eletyp, elecon)
    if u is None:
        f.write_vtu_file(**kwds)
    else:
        kw = dict([(key, zeros_like(val)) for (key, val) in kwds.items()])
        f.snapshot(**kw)
        f.snapshot(u=u, **kwds)

def ReadMesh(filename, disp=0):
    f = VTKFile(filename=filename, mode='r')
    if disp:
        return f.coord, f.nodlab, f.elecon, f.elelab

    # Create node and element tables
    nodtab = [[f.nodlab[i]]+xc.tolist() for i,xc in enumerate(f.coord)]

    numele = len(f.elecon)
    eletab = [[f.elelab[e]]+f.elecon[e].tolist() for e in range(numele)]

    return nodtab, eletab, f.nodesets, f.elemsets, f.surfaces

def test_write_fe_results():
    from distmesh import drectangle, distmesh2d, huniform
    random.seed(190) # Always the same results
    fd = lambda p: drectangle(p, -1, 1, -1, 1)
    fh = huniform
    coord, elecon = distmesh2d(fd, fh, .1, (-1, -1, 1, 1),
                               [(-1, -1),(-1, 1),(1, -1),(1, 1)])
    jobid = 'Job'
    nodlab = range(coord.shape[0])
    nodmap = dict([(n,n) for n in nodlab])
    elelab = range(elecon.shape[0])
    elemap = dict([(n,n) for n in elelab])
    eletyp = [ElementFamily(2,3)] * elecon.shape[0]
    scal = random.rand(coord.shape[0])
    vect = random.rand(coord.shape[0]*2).reshape(-1,2)
    tens = random.rand(elecon.shape[0]*9).reshape(-1,9)
    symt = random.rand(elecon.shape[0]*6).reshape(-1,6)
    kwds = dict(scal=scal,vect=vect,tens=tens,symt=symt)
    u = zeros_like(coord)
    u[:,0] = 1
    WriteFEResults(jobid, coord, nodmap, elemap, eletyp, elecon, u=u, **kwds)
    filename = jobid+'.vtu'
    WriteVTUMesh(filename, coord, nodlab, elelab, eletyp, elecon, check=1)
    os.remove(filename)

if __name__ == '__main__':
    #ReadMesh('uniform_plate_tri_0.05.vtu')
    test_write_fe_results()
