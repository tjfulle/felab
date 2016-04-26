import os
from numpy import *
import logging
from argparse import ArgumentParser
from collections import OrderedDict

from ..constants import *
from ..utilities import is_listlike, UserInputError
from ..elemlib import ElementFamily
from ..material import Material

__all__ = ['Mesh', 'UnitSquareMesh', 'RectilinearMesh2D', 'VTKMesh',
           'AbaqusMesh', 'GenesisMesh']

def is_listlike(a):
    return (not hasattr(a, 'strip') and
            hasattr(a, '__getitem__') or hasattr(a, '__iter__'))

def is_stringlike(s):
    return hasattr(s, 'strip')

class Container(OrderedDict):
    def __contains__(self, key):
        return key.lower() in [k.lower() for k in self.keys()]
    def get(self, key, default=None):
        K = key.lower()
        for (k, v) in self.items():
            if k.lower() == K:
                return v
        return default

class ElementBlock:
    def __init__(self, name, id, labels, elefam, elecon):
        self.name = name.upper()
        self.id = id
        self.labels = labels
        self.elefam = elefam
        self.numele = len(labels)
        self.elecon = elecon
        self.material = None
    @property
    def eletyp(self):
        return self.elefam
    @eletyp.setter
    def eletyp(self, arg):
        self.elefam = arg
    def set_material(self, material):
        if self.material is not None:
            if material == self.material:
                return None
            logging.warn('CHANGING MATERIAL MODEL')
        if not isinstance(material, Material):
            raise UserInputError('NOT A MATERIAL')
        self.material = material

class Mesh(object):
    """Class for creating a finite element mesh.

    Parameters
    ----------
    nodtab : list of list
        ``nodtab[n]`` is a list describing the nth node where
        ``nodtab[n][0]`` is the node label of the nth node and
        ``nodtab[n][1:]`` are the nodal coordinates of the nth node
    eletab : list of list
        ``eletab[e]`` is a list describing the eth element where
        ``eletab[e][0]`` is the element label of the eth element and
        ``eletab[e][1:]`` are the nodes forming the eth element, the nodal
        labels must correpond to the node labels in ``nodtab``
    filename : str
        Name of a genesis or vtu file containing the mesh
    p : ndarray of float
        Array of nodal coordinates
    t : ndarray of int
        Array of element connectivity

    Notes
    -----
    The nodtab/eletab, filename, and p/t arguments are mutually exclusive.

    If ``p`` and ``t`` are specified, node labels are assigned to nodes in the
    order they appear in the nodal coordinates array ``p`` and take values of
    ``0`` to ``n-1``, where ``n`` is the length of ``p``. The element
    connectivity array ``t`` must use this numbering scheme. This method is
    useful for creating meshes using triangulations created by the
    ``distmesh2d`` module.

    Examples
    --------

    In the following examples, a simple mesh is created of a unit square
    centered at the origin with node and element numbers as follows:

    .. image:: unitsquare1.png
       :scale: 35
       :align: center

    Method 1: specifying ``nodtab`` and ``eletab``:

    >>> nodtab = [[10, -.5, -.5], [20, .5, -.5], [30, .5, .5], [40, -.5, .5]]
    >>> eletab = [[100, 10, 20, 30, 40]]
    >>> mesh = Mesh(nodtab=nodtab, eletab=eletab)

    Method 2: specifying ``filename``

    >>> mesh = Mesh(filename='unitsquare.g')

    Method 3: specifying ``p`` and ``t``:

    >>> p = array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
    >>> t = array([[0, 1, 2, 3]])
    >>> mesh = Mesh(p=p, t=t)

    """
    type = 'default'
    def __init__(self, nodtab=None, eletab=None, filename=None, p=None, t=None):
        self._bndry_nod2 = None
        self._free_edges = None
        self.element_blocks = Container()

        o1 = nodtab is not None and eletab is not None
        o2 = filename is not None
        o3 = p is not None and t is not None
        if not (o1 or o2 or o3):
            raise UserInputError('nodtab/eletab, filename, or p/t required')
        elif len([x for x in (o1, o2, o3) if x]) > 1:
            raise UserInputError('nodtab/eletab, filename, an p/t are '
                                 'mutually exclusive')
        if o1:
            # Format nodes and elements
            data = self.parse_nod_and_elem_tables(nodtab, eletab)
            self.init1(*data)

        elif o2:
            self.init_from_file(filename)

        elif o3:
            p, t = asarray(p, dtype=float), asarray(t, dtype=int)
            nodmap = dict(zip(range(p.shape[0]), range(p.shape[0])))
            eletab = dict([(i,t[i]) for i in range(t.shape[0])])
            self.init1(nodmap, p, eletab)

        else:
            raise UserInputError('nodtab/eletab or filename required')

    @property
    def unassigned(self):
        return self.num_assigned != self.numele

    def init1(self, nodmap, coord, eletab, nodesets=None,
              elemsets=None, surfaces=None):
        self.nodmap = nodmap
        self.coord = asarray(coord)
        self.numnod, self.dimensions = self.coord.shape

        self.eletab = eletab
        self.numele = len(self.eletab)
        self.elemap = dict(zip(self.eletab.keys(), [None]*self.numele))
        self.ielemap = {}

        self.nodesets = nodesets or {}
        self.elemsets = elemsets or {}
        self.surfaces = surfaces or {}

        self.eleblx = []

        # maximum number of edges on any one element
        self.maxedge = max([len(ElementFamily(self.dimensions,len(nodes)).edges)
                            for nodes in self.eletab.values()])

        # Number of elements assigned to a block
        self.num_assigned = 0

    def init2(self, nodmap, coord, elemap, eleblx, nodesets, elemsets, sidesets):
        self.nodmap = nodmap
        self.coord = asarray(coord)
        self.numnod, self.dimensions = self.coord.shape

        self.elemap = elemap
        self.numele = len(self.elemap)
        self.ielemap = array(sorted(elemap.keys(), key=lambda k: elemap[k]))

        self.nodesets, self.elemsets, self.surfaces = {}, {}, {}
        for (name, nodes) in nodesets.items():
            self.nodesets[name.upper()] = nodes
        for (name, elems) in elemsets.items():
            self.elemsets[name.upper()] = elems
        for (name, surf) in sidesets.items():
            self.surfaces[name.upper()] = surf

        self.eleblx = eleblx
        for eb in eleblx:
            self.element_blocks[eb.name.upper()] = eb
        self.maxedge = max([len(eb.elefam.edges) for eb in eleblx])
        self.num_assigned = self.numele

    def init_from_file(self, filename):
        import pyfem2.mesh.aba as aba
        import pyfem2.mesh.vtk as vtk
        import pyfem2.mesh.exodusii as exodusii
        if filename.endswith(('.vtk', '.vtu')):
            data = self.parse_nod_and_elem_tables(*vtk.ReadMesh(filename))
            self.init1(*data)

        elif filename.endswith(('.exo', '.g', '.gen')):
            exo = exodusii.File(filename, mode='r')
            self.init2(exo.nodmap, exo.coord, exo.elemap, exo.eleblx,
                       exo.nodesets, exo.elemsets, exo.sidesets)

        elif filename.endswith('.inp'):
            self.type = 'abaqus'
            data = aba.ReadInput(filename)
            nodtab, eletab, nodesets, elemsets, surfaces, eleblx = data
            nodmap, coord, eletab1 = self.parse_nod_and_elem_tables(nodtab, eletab)
            self.init1(nodmap, coord, eletab1)
            for (name, (et, elems, ebmat)) in eleblx.items():
                self.ElementBlock(name, elems, material=ebmat)
            for (name, nodes) in nodesets.items():
                self.NodeSet(name, nodes)
            for (name, elems) in elemsets.items():
                self.ElementSet(name, elems)
            for (name, surf) in surfaces.items():
                self.Surface(name, surf)

        else:
            raise UserInputError('Unknown file type')

    def get_internal_node_ids(self, label):
        if is_stringlike(label) and label.upper() in self.nodesets:
            return self.nodesets[label.upper()]
        elif isinstance(label, int):
            inodes = [self.nodmap[label]]
        elif label == ALL:
            inodes = range(self.numnod)
        elif label == BOUNDARY:
            inodes = self.boundary_nodes()
        elif label in (ILO, IHI, JLO, JHI, KLO, KHI):
            inodes = self.nodes_in_rectilinear_region(label)
        elif is_stringlike(label):
            raise UserInputError('No such node set {0!r}'.format(label))
        else:
            inodes = [self.nodmap[xn] for xn in label]
        return array(inodes, dtype=int)

    @staticmethod
    def parse_nod_and_elem_tables(nodtab, eletab):
        """

        .. _parse_nodes_and_elements:

        Format the node map, coordinates array, element map, types, and
        connectivity table

        Parameters
        ----------
        nodtab : list of list
            ``nodtab[n]`` is a list describing the nth node where
            ``nodtab[n][0]`` is the node label of the nth node and
            ``nodtab[n][1:]`` are the nodal coordinates of the nth node
        eletab : list of list
            ``eletab[e]`` is a list describing the eth element where
            ``eletab[e][0]`` is the element label of the eth element and
            ``eletab[e][1:]`` are the nodes forming the eth element, the
            nodal labels must correpond to the node labels in ``nodtab``

        Returns
        -------
        nodmap : dict
            nodmap[xn] is the internal node number of node labeled xn
        coord : ndarray
            Nodal coordinates.  coord[:,0] are the x, coord[:,1] the y, etc.
        eletab1 : dict
            eletab1[xel] are the internal node numbers of element labeled xel

        """
        # Basic node info
        nodmap = {}
        numnod, maxdim = len(nodtab), 0
        for (inode, noddef) in enumerate(nodtab):
            if noddef[0] in nodmap:
                raise UserInputError('DUPLICATE NODE LABEL: {0}'.format(noddef[0]))
            nodmap[noddef[0]] = inode
        maxdim = max(maxdim, len(noddef[1:]))
        coord = zeros((numnod, maxdim))
        for (inode, noddef) in enumerate(nodtab):
            coord[inode,:len(noddef[1:])] = noddef[1:]

        # Basic element info
        eletab1 = {}
        for (iel, eledef) in enumerate(eletab):
            if eledef[0] in eletab1:
                raise UserInputError('Duplicate element label: '
                                     '{0}'.format(eledef[0]))
            eletab1[eledef[0]] = [nodmap[n] for n in eledef[1:]]

        return nodmap, coord, eletab1

    def find_edges(self):
        # Find edges
        if self.unassigned:
            raise UserInputError('Mesh element operations require all elements be '
                                 'assigned to an element block')
        if self._free_edges is not None:
            return self._free_edges

        edges = {}
        k = 0
        for (ieb, eb) in enumerate(self.eleblx):
            for (j, elenod) in enumerate(eb.elecon):
                iel = self.elemap[eb.labels[j]]
                for (i, ix) in enumerate(eb.elefam.edges):
                    edge = tuple(elenod[ix])
                    edges.setdefault(tuple(sorted(edge)), []).append((iel, i)+edge)
        self._free_edges = edges.values()
        return self._free_edges

    def boundary_nodes2(self):
        if self._bndry_nod2 is not None:
            return self._bndry_nod2
        edges = self.find_edges()
        nodes = []
        for edge in edges:
            if len(edge) == 1:
                # edge not shared by multiple elements
                nodes.extend(edge[0][2:])
        self._bndry_nod2 = unique(sorted(nodes))
        return self._bndry_nod2

    def boundary_nodes(self):
        if self.dimensions == 1:
            return [argmin(self.coord), argmax(self.coord)]
        elif self.dimensions == 2:
            return self.boundary_nodes2()
        raise UserInputError('3D meshes not supported')

    def nodes_in_rectilinear_region(self, region, tol=1e-6):
        # Find nodes in region
        if region not in (ILO, IHI, JLO, JHI, KLO, KHI):
            raise UserInputError('unknown region {0!r}'.format(region))
        axis, fun = {ILO: (0, amin),
                     IHI: (0, amax),
                     JLO: (1, amin),
                     JHI: (1, amax),
                     KLO: (2, amin),
                     KHI: (2, amax)}[region]
        xpos = fun(self.coord[:,axis])
        return where(abs(self.coord[:,axis] - xpos) < tol)[0]

    def find_surface(self, region):
        if self.unassigned:
            raise UserInputError('Mesh element operations require all elements be '
                                 'assigned to an element block')
        if is_stringlike(region) and region.upper() in self.surfaces:
            return self.surfaces[region.upper()]
        if region in (ILO, IHI, JLO, JHI, KLO, KHI):
            if self.dimensions == 1:
                return self.find_surface1(region)
            if self.dimensions == 2:
                return self.find_surface2(region)
        elif is_stringlike(region):
            raise UserInputError('No such surface {0!r}'.format(region))
        # Check if region is a surface
        if not is_listlike(region):
            raise UserInputError('Unrecognized surface: {0}'.format(region))
        surface = []
        region = asarray(region)
        if region.ndim == 1:
            region = region[newaxis,:]
        for item in region:
            elem, edge = [int(x) for x in item]
            if is_listlike(elem):
                surface.extend([(self.elemap[e], edge) for e in elem])
            else:
                surface.append((self.elemap[elem], edge))
        return surface

    def find_surface1(self, region):
        if region not in (ILO, IHI):
            raise UserInputError('INCORRECT 1D REGION')
        if region == ILO:
            node = argmin(self.coord)
        else:
            node = argmax(self.coord)
        for (ieb, eb) in enumerate(self.eleblx):
            for (e, elenod) in enumerate(eb.elecon):
                if node in elenod:
                    i = 0 if node == elenod[0] else 1
                    return [self.elemap[eb.labels[e]], i]

    def find_surface2(self, region):
        nodes = self.nodes_in_rectilinear_region(region)
        surface = []
        for (ieb, eb) in enumerate(self.eleblx):
            for (e, elenod) in enumerate(eb.elecon):
                w = where(in1d(elenod, nodes))[0]
                if len(w) < 2:
                    continue
                w = tuple(sorted(w))
                for (ie, edge) in enumerate(eb.elefam.edges):
                    if tuple(sorted(edge)) == w:
                        # the internal node numbers match, this is the edge
                        surface.append((self.elemap[eb.labels[e]], ie))
        return array(surface)

        for (e, c) in enumerate(self.elecon):
            c = c[:self.elefam[e].nodes]
            w = where(in1d(c, nodes))[0]
            if len(w) < 2:
                continue
            w = tuple(sorted(w))
            for (ie, edge) in enumerate(self.elefam[e].edges):
                if tuple(sorted(edge)) == w:
                    # the internal node numbers match, this is the edge
                    surface.append((e, ie))
        return array(surface)

    def NodeSet(self, name, region):
        if not is_stringlike(name):
            raise UserInputError('NAME MUST BE A STRING')
        self.nodesets[name.upper()] = self.get_internal_node_ids(region)

    def Surface(self, name, surface):
        if self.unassigned:
            raise UserInputError('MESH ELEMENT OPERATIONS REQUIRE ALL ELEMENTS BE '
                                 'ASSIGNED TO AN ELEMENT BLOCK')
        if not is_stringlike(name):
            raise UserInputError('NAME MUST BE A STRING')
        self.surfaces[name.upper()] = self.find_surface(surface)

    def ElementSet(self, name, region):
        if self.unassigned:
            raise UserInputError('MESH ELEMENT OPERATIONS REQUIRE ALL ELEMENTS BE '
                                 'ASSIGNED TO AN ELEMENT BLOCK')
        if not is_stringlike(name):
            raise UserInputError('NAME MUST BE A STRING')
        if region == ALL:
            ielems = range(self.numele)
        else:
            if not is_listlike(region):
                region = [region]
            ielems = [self.elemap[el] for el in region]
        self.elemsets[name.upper()] = array(ielems, dtype=int)

    def ElementBlock(self, name, elements, material=None):

        if name in self.element_blocks:
            raise UserInputError('{0!r} ALREADY AN ELEMENT BLOCK'.format(name))
        if elements == ALL:
            xelems = sorted(self.elemap.keys())
        else:
            if not is_listlike(elements):
                elements = [elements]
            xelems = elements

        if not self.unassigned:
            raise UserInputError('ALL ELEMENTS HAVE BEEN ASSIGNED')

        # Elements are numbered in the order the appear in element blocks,
        # we need to map from the temporary internal numbers
        blkcon, badel = [], []
        numblkel = len(xelems)
        ielems = arange(self.num_assigned, self.num_assigned+numblkel)
        for (i, xelem) in enumerate(xelems):
            # external element label -> new internal element ID
            ielem = ielems[i]
            self.elemap[xelem] = ielem
            self.ielemap[ielem] = xelem
            elenod = self.eletab[xelem]
            if blkcon and len(elenod) != len(blkcon[-1]):
                badel.append(xelem)
                continue
            blkcon.append(elenod)
        if badel:
            badel = ',\n   '.join(str(e) for e in badel)
            raise UserInputError('THE FOLLOWING ELEMENTS HAVE INCONSISTENT ELEMENT '
                                 'CONNECTIVITY:\n   {0}'.format(badel))
        blkcon = array(blkcon, dtype=int)
        elefam = ElementFamily(self.dimensions, blkcon.shape[1])
        blk = ElementBlock(name, len(self.eleblx)+1, xelems, elefam, blkcon)
        self.eleblx.append(blk)
        self.element_blocks[blk.name] = blk
        self.num_assigned += len(ielems)
        self.maxedge = max(self.maxedge, len(elefam.edges))

        if material is not None:
            blk.set_material(material)
        return blk

    def to_genesis(self, filename):
        exof = File(filename, mode='w')
        if not self.eleblx:
            # put in a single element block
            self.ElementBlock('ElementBlock1', ALL)
        exof.genesis(self.nodmap, self.elemap, self.coord,
                     self.eleblx, nodesets=self.nodesets, elemsets=self.elemsets,
                     sidesets=self.surfaces)
        exof.close()

    def Plot2D(self, xy=None, elecon=None, u=None, color=None, ax=None, show=0,
               weight=None, colorby=None, linestyle='-', label=None, xlim=None,
               ylim=None, filename=None, **kwds):
        assert self.dimensions == 2
        from matplotlib.patches import Polygon
        import matplotlib.lines as mlines
        from matplotlib.collections import PatchCollection
        from matplotlib.cm import coolwarm, Spectral
        import matplotlib.pyplot as plt

        if xy is None:
            xy = array(self.coord)
        if elecon is None:
            elecon = []
            for blk in self.eleblx:
                elecon.extend(blk.elecon.tolist())
            elecon = asarray(elecon)
        if u is not None:
            xy += u.reshape(xy.shape)

        patches = []
        for points in xy[elecon[:]]:
            quad = Polygon(points, True)
            patches.append(quad)

        if ax is None:
            fig, ax = plt.subplots()

        #colors = 100 * random.rand(len(patches))
        p = PatchCollection(patches, linewidth=weight, **kwds)
        if colorby is not None:
            colorby = asarray(colorby).flatten()
            if len(colorby) == len(xy):
                # average value in element
                colorby = array([average(colorby[points]) for points in elecon])
            p.set_cmap(Spectral)  #coolwarm)
            p.set_array(colorby)
            p.set_clim(vmin=colorby.min(), vmax=colorby.max())
            fig.colorbar(p)
        else:
            p.set_edgecolor(color)
            p.set_facecolor('None')
            p.set_linewidth(weight)
            p.set_linestyle(linestyle)

        if label:
            ax.plot([], [], color=color, linestyle=linestyle, label=label)

        ax.add_collection(p)

        if not ylim:
            ymin, ymax = amin(xy[:,1]), amax(xy[:,1])
            dy = max(abs(ymin*.05), abs(ymax*.05))
            ax.set_ylim([ymin-dy, ymax+dy])
        else:
            ax.set_ylim(ylim)

        if not xlim:
            xmin, xmax = amin(xy[:,0]), amax(xy[:,0])
            dx = max(abs(xmin*.05), abs(xmax*.05))
            ax.set_xlim([xmin-dx, xmax+dx])
        else:
            ax.set_xlim(xlim)
        ax.set_aspect('equal')

        if show:
            plt.show()

        if filename is not None:
            plt.legend()
            plt.savefig(filename, transparent=True,
                        bbox_inches="tight", pad_inches=0)

        return ax

    def PlotScalar2D(self, u, show=0):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib.cm import Spectral
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        elecon = []
        for eb in self.eleblx:
            elecon.extend(eb.elecon)
        elecon = asarray(elecon)
        ax.plot_trisurf(self.coord[:,0], self.coord[:,1], u,
                        triangles=elecon, cmap=Spectral)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        if show:
            plt.show()
        return

    def PutNodalSolution(self, filename, u):
        import pyfem2.mesh.exodusii as exodusii
        if not self.eleblx:
            self.ElementBlock('ElementBlock1', ALL)
        if not filename.endswith(('.exo', '.e')):
            filename += '.exo'
        exodusii.PutNodalSolution(filename, self.nodmap, self.elemap, self.coord,
                                  self.eleblx, u)

def GenesisMesh(filename):
    """
    Generates a finite element mesh from a Genesis file.

    Parameters
    ----------
    filename : str
        The path to a valid Genesis file

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls ``mesh.Mesh`` with the ``filename`` keyword and
    stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

    """
    if not os.path.isfile(filename):
        raise UserInputError('NO SUCH FILE {0!r}'.format(filename))
    return Mesh(filename=filename)

def AbaqusMesh(filename):
    """
    Generates a finite element mesh from a Abaqus input file.

    Parameters
    ----------
    filename : str
        The path to a valid Genesis file

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls ``mesh.Mesh`` with the ``filename`` keyword and
    stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

    """
    if not os.path.isfile(filename):
        raise UserInputError('NO SUCH FILE {0!r}'.format(filename))
    return Mesh(filename=filename)

def VTKMesh(filename):
    """
    Generates a finite element mesh from a vtk .vtu file.

    Parameters
    ----------
    filename : str
        The path to a valid .vtu file

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls ``mesh.Mesh`` with the ``filename`` keyword and
    stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

    """
    if not os.path.isfile(filename):
        raise UserInputError('NO SUCH FILE {0!r}'.format(filename))
    return Mesh(filename=filename)

def RectilinearMesh2D(nx=1, ny=1, lx=1, ly=1, shift=None):
    """
    Generates a rectilinear 2D finite element mesh.

    Parameters
    ----------
    shape : tuple
        (nx, ny) where nx is the number elements in :math:`x` and ny
        number of element in :math:`y`.
    lengths : tuple
        (lx, ly) where lx is the length of the mesh in :math:`x` and ny
        is the length of the mesh in :math:`y`.

    Returns
    -------
    Mesh object

    """
    nodtab, eletab = GenRectilinearMesh2D((nx, ny), (lx, ly), shift=shift)
    return Mesh(nodtab=nodtab, eletab=eletab)

def UnitSquareMesh(nx=1, ny=1, shift=None):
    """
    Generates a rectilinear 2D finite element mesh.

    Parameters
    ----------
    shape : tuple
        (nx, ny) where nx is the number elements in :math:`x` and ny
        number of element in :math:`y`.

    Returns
    -------
    Mesh object

    Notes
    -----
    This method calls the ``Mesh.RectilinearMesh2D`` class method and
    stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

    """
    nodtab, eletab = GenRectilinearMesh2D((nx, ny), (1, 1), shift=shift)
    return Mesh(nodtab=nodtab, eletab=eletab)

def GenRectilinearMesh2D(shape, lengths, shift=None):
    nx, ny = shape
    if nx < 1:
        raise UserInputError('REQURES AT LEAST 1 ELEMENT IN X')
    if ny < 1:
        raise UserInputError('REQURES AT LEAST 1 ELEMENT IN Y')

    if shift is None:
        shift = zeros(2)
    else:
        shift = asarray(shift)

    shape = asarray([nx+1,ny+1])
    lx, ly = lengths
    xpoints = linspace(0, lx, nx+1) + shift[0]
    ypoints = linspace(0, ly, ny+1) + shift[1]
    coord = array([(x, y) for y in ypoints for x in xpoints])
    numnod = prod(shape)
    numele = prod(shape - 1)

    # Connectivity
    row = 0
    elecon = zeros((numele, 4), dtype=int)
    nelx = xpoints.size - 1
    for elem_num in range(numele):
        ii = elem_num + row
        elem_nodes = [ii, ii + 1, ii + nelx + 2, ii + nelx + 1]
        elecon[elem_num, :] = elem_nodes
        if (elem_num + 1) % (nelx) == 0:
            row += 1
        continue

    nodtab = []
    for n in range(numnod):
        nodtab.append([n+1])
        nodtab[-1].extend(coord[n])
    eletab = []
    for e in range(numele):
        eletab.append([e+1])
        eletab[-1].extend([n+1 for n in elecon[e]])

    return nodtab, eletab

def VTU2Genesis(nodtab=None, eletab=None, filename=None):
    if filename is None:
        assert nodtab is not None and eletab is not None
        outfile = 'mesh.g'
    elif not os.path.isfile(filename):
        assert nodtab is not None and eletab is not None
        assert filename.endswith('.g')
        outfile = filename
        filename = None
    else:
        assert nodtab is None and eletab is None
        outfile = os.path.splitext(filename)[0] + '.g'
    try:
        mesh = Mesh(nodtab=nodtab, eletab=eletab, filename=filename)
    except KeyError:
        return
    mesh.to_genesis(outfile)

def INP2Genesis(filename):
    lines = open(filename).readlines()
    kw, name = None, None
    nodtab = []
    eletab = []
    nodesets = {}
    elemsets = {}
    eleblx = {}
    for line in lines:
        line = ','.join([x.strip() for x in line.split(',')])
        if line.startswith('**'):
            continue
        if not line.split():
            continue
        if line.startswith('*'):
            name = None
            line = line.split(',')
            kw = line[0][1:].upper()
            opts = {}
            for opt in line[1:]:
                k, v = opt.split('=')
                opts[k.strip().upper()] = v.strip()
            if kw != 'ELEMENT':
                name = None
            if kw == 'ELEMENT':
                name = opts.get('ELSET')
                if name is None:
                    raise UserInputError('REQUIRES ELEMENTS BE PUT IN ELSET')
                eleblx[name.upper()] = []
            elif kw == 'NSET':
                name = opts['NSET']
                nodesets[name.upper()] = []
            elif kw == 'ELSET':
                elemsets[name.upper()] = []
            continue
        if kw is None:
            continue
        if kw == 'NODE':
            line = line.split(',')
            nodtab.append([int(line[0])] + [float(x) for x in line[1:]])
            continue
        elif kw == 'ELEMENT':
            eledef = [int(n) for n in line.split(',')]
            eletab.append(eledef)
            eleblx[name].append(eledef[0])
            continue
        elif kw == 'NSET':
            nodesets[name.upper()].extend([int(n) for n in line.split(',')
                                           if n.split()])
            continue
        elif kw == 'ELSET':
            elemsets[name.upper()].extend([int(n) for n in line.split(',')
                                           if n.split()])
            continue
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    for (name, elems) in eleblx.items():
        mesh.ElementBlock(name, elems)
    for (name, nodes) in nodesets.items():
        mesh.NodeSet(name, nodes)
    for (name, elems) in elemsets.items():
        mesh.ElementSet(name, elems)
    outfile = os.path.splitext(filename)[0] + '.g'
    mesh.to_genesis(outfile)
    return
